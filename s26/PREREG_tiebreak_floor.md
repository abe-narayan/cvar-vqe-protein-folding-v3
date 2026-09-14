# PREREG_tiebreak_floor -- THE PRODUCTION ENDPOINT'S NOISE FLOOR FROM BLOSUM TIE-BREAKING AT THE K = 500 BOUNDARY (lane W, Sprint 26)

Written 2026-09-13 09:30, before any result. Idea: `s26/IDEA_tiebreak_noise_floor.md`. Code:
`s26/w_tiebreak.py` (synthetic tests `s26/w_tiebreak_test.py`). Results:
`s26/results/w_tiebreak_*.json`. Never edited after the first gated run; addenda appended. The
native-free half (`draws`) reads no native; the gated half (`endpoint`) refuses to run until
"PHASE 0 SIGNED OFF" is in `s26/LEDGER.md`. Runs only if the tournament ranks it or the
coordinator asks; the one-target probe runs now so it can start at once.

## 0. What the record knows

- The pool is the top 500 windows by BLOSUM62 sum (an integer) under a stable argsort over the
  universe's corpus order; the `pdbs/` file set is pinned because "sort order feeds BLOSUM
  tie-breaking" (state brief section 3).
- Native-free census (`s26/results/w_selfcopy_census.json`, L30): 126/126 targets have a tie
  at the boundary; median 115 windows share the boundary score, about 57 inside and 60 outside
  the pool; so about 11% of every pool is set by corpus order.
- S17 L628: the ORACLE pool best under random tie-breaks is 1.708 A against 1.711 (sd 0.018);
  the 11.4% figure was called robust for the oracle best. S8-12 / S10 (`docs/FINDINGS.md:3668,
  6106`): "a different tie-break gives a different pool" (1A1P pool best 3.270 vs 3.336). S9
  (`docs/FINDINGS.md:1526, 7261`): 0.004 to 0.005 A residuals on 126-target means between the
  stable and the plain argsort, attributed to tie-break noise. Nobody measured the built-chain
  endpoint's spread over random tie-breaks or the paired MDE between two of them.
- Recorded effects at the hundredths level that a floor could sit under: restrained AMBER
  relaxation of the average -0.022 [-0.036, -0.009] (memory `averaging-space-beats-the-objective`),
  the ORACLE functional weight 0.015 (S24 L16), all-atom reranking +0.004 (memory
  `nothing-ranks-within-the-pool`), C27's +0.0004, the projection-mode differences of 1e-2 (L19).

## 1. Hypothesis and exact falsifier

H_T: re-drawing the boundary members of the pool at random (a uniformly random subset of the
boundary tie class of the same size as the production convention's, the non-tied part of the
pool unchanged, everything downstream identical) moves the built chain by a per-target sd s_tie
and the 126-target mean by m_tie, and the paired MDE between two random draws is the smallest
effect the pipeline can distinguish from its own unstated convention.

Falsifier. The claim of the idea ("at least one recorded effect is inside the floor") is
FALSIFIED if m_tie < 0.002 A AND the paired MDE between two draws < 0.010 A on the built chain.
It is CONFIRMED if any of the recorded effects listed in section 0 is below the paired MDE
between two draws. Registered predictions: m_tie 0.003 to 0.010 A; median s_tie 0.03 to 0.08
A; paired MDE 0.02 to 0.05 A; about 8 of the 75 averaged members change per draw. Secondary,
reported not decided on: (i) whether the production convention (corpus order) differs from
the mean of the draws beyond that MDE (a sign would be a finding about the pinned `pdbs/`
order); (ii) the same floor on the `sel` basis (argmin over the pool) and the point cloud.

## 2. Design

Per target: the boundary score s = sim[order[499]]; the tie class T = {w : sim[w] == s};
n_in = |T inside the production pool|. Draw k (k = 0..7): `s15.seed.stable_rng("w_tiebreak",
pdb, k)` chooses n_in members of T uniformly without replacement; the pool is the non-tied
members (sim > s, unchanged) followed by the drawn members in corpus order. Then the shipped
score, top-75, medoid-frame average and `I.project` (ramah 0.3, multi-start), through
`s26/w_selfcopy.emit`, which reproduces the production top-75 and cloud (gate: `sub` equal,
`avg_ca` to 9e-16, chain 1.1e-4 A after superposition). Native-free outputs per draw: the
emissions (cloud, chain, fit), the fraction of the top-75 replaced, the argmin identity, and the
triangle bounds against the production emission (RMSD between chains bounds the change in
RMSD-to-native). Gated outputs: RMSD-to-native per draw and basis; s_tie, m_tie; `ST.compare`
between draw 0 and draw 1, and between the production emission and the mean over draws, with
fold-clustered CIs; the distribution of paired MDEs over the 28 draw pairs.

Controls. The zero-change control is draw = production convention (the gate). The matched
control for "any change of 8 members" is not needed: the draws are that operator. Ties are
handled by the tie set (`argmin_set`) and `ST.argmin_tied`'s rule for `sel`.

## 3. Expected effect against the MDE

This experiment measures an MDE. Its numbers are the deliverable.

## 4. Memory, time, agent-hours

`s26/w_selfcopy.py`'s retrieval probe ran at 0.104 GB for one target with three projections;
one universe at a time. 126 targets x 8 draws x ~3 s = 50 min CPU, checkpointed per target,
est-ram 0.5 GB, tag CPU. Probe: one target, 8 draws, under jobrun, before the full run.
Agent-hours: 1 code, 1 run, 1 write-up.

## 5. Operator forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| which ties | the boundary class only (the only class that changes membership) | interior ties (they change the pool's ORDER, not its membership; the shipped score is order-blind except for float ties in the top-75 cut) |
| draw law | uniform subsets of the tie class, same size as production | weighting by anything |
| number of draws | 8 (28 paired MDEs) | 16 (doubles the cost; add as an addendum if the 8-draw sd of m_tie is itself uncertain by more than 30%) |
| basis | built chain PRIMARY; cloud and sel carried | AMBER-relaxed emission |

## ADDENDUM 1 (2026-09-13 19:45) -- the one-target probe; nothing above edited; no launch

- Synthetic tests `s26/w_tiebreak_test.py`: ALL OK (job `w_tiebreak_test`, 5 s).
- Probe `w_tiebreak_probe` (1CEK, 8 draws; exit 0, 35 s, peak RSS 0.098 GB;
  `s26/results/w_selfcopy_tiebreak_probe_1CEK.json`): boundary tie class 130 windows, 44 inside
  the production pool; production gate true (top-75 == `sub`); per draw the pool overlap with the
  production pool is 0.932 to 0.946, 3 to 7 of the 75 averaged members change, the argmin never
  changes, and the built chain moves 0.008 / 0.030 / 0.024 / 0.023 / 0.021 / 0.011 / 0.031 / 0.008 A
  (triangle bounds; the lam = 0 chain the same to 1e-3). The registered prediction "about 8 of 75
  change" is 3 to 7 on this target.
- Scale of the full run from the probe: 126 x 8 x ~3.5 s = 50 min CPU, est-ram 0.5 GB.
- The L44 control population (one member of the 75 replaced) already shows the two discrete
  operators this floor would price: the medoid frame and the projection branch, with signed
  moves up to 0.51 A on the built chain on 5% of targets.
- Launch condition: only if the Adversary ranks the idea and the coordinator confirms headroom
  (L42: one sub-1 GB job at a time). Not launched.
