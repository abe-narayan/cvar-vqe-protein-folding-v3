# PREREG_S29_D_m6 -- THE FIXED-PROFILE CONTROL: REPLACE THE ENTIRE QUANTUM STAGE BY A TARGET-INDEPENDENT RANK-WEIGHT PROFILE

Lane D (Adversary), Sprint 29. Written 2026-09-20 00:44 Pacific, **before any number of this
experiment has been read**. Commissioned by the coordinator from lane T's S29-L15 clause (M6).
Contract `s29/S29_CONTRACT.md` rules 5, 6, 7, 9, 11, 15, 16, 17, 18. Code `s29/s29_D_m6.py`;
results `s29/results/s29_D_m6*.json`; tests `tests/test_s29_D.py`.

**Declaration of what I had already seen when this was written** (it is the only thing that can
compromise this pre-registration, and it is why the falsifier below is stated in the direction it
is): the schema and summary statistics of `s27/results/vqe_rows.jsonl`, namely that the DIS rows
exist for 126 targets at seeds 0 and 1, that the realised tail size m has mean 74.07 / sd 6.77
(seed 0) and 71.19 / 7.97 (seed 1), min 52 max 85 / min 49 max 82, that m equals exactly 74 on 3
targets and exactly 75 on 9 (seed 0), and that each row carries `rmsd_vqe`, `rmsd_topm` and
`rmsd_top75`. I have NOT computed any contrast, any per-target difference, or any chain number.

## 0. Labels and stakes
Nothing here is ORACLE-tuned; the RMSDs are the endpoint (ORACLE by definition, as every RMSD is)
and no native quantity chooses any parameter. This is the charter's central requirement --
"removing or randomising the quantum stage measurably degrades the result" -- tested at its
sharpest on the DEPLOYED spine, and the S29 contract already concedes the shape of the answer
(the production spine "is exactly that" system). The experiment's value is to state it with
numbers, on the reporting basis, with the honest framing lane T supplied: this control does not
diminish the trainability result (the optimiser does reduce the objective, S28-L21, S29-L15); it
measures what that reduction buys at the endpoint.

## 1. The control, defined so that it contains no per-target computation
Lane T's derivation (S29-L15): the deployed CVaR is exactly constant along D - m - 1 = 437 of 511
simplex directions; the entropy term makes the optimum exactly uniform above the VaR and
exponentially enhanced below it; the energies are `zrank` of the score, the same standardised
rank ladder on every target to 1.18% of range (S25 L17); therefore the optimal p* depends on
(alpha, T) and nothing else, and the emitted structure's only target-specific input is which
candidate the distogram put at which rank.
**FIXED-PROFILE CONTROL.** Compute p*(alpha = 0.18, T = 0.5) ONCE on the standardised rank ladder
(no target, no circuit, no optimiser). Apply it to each target's own score-sorted candidate list.
Read it out with the deployed readout (`s24.d_harness.readout_uniform` over the CVaR tail set),
which by the set-equality theorem is the uniform average of the top-m* prefix, m* being the
prefix at which p*'s cumulative mass reaches alpha. m* is one number for all 126 targets.
**Arms, all registered now** (m* from p*, plus a sweep so the answer is a curve and not one cell):
m* in {the computed p* prefix, 70, 71, 74, 75, 80} -- 75 is production itself, which makes the
production anchor a member of the control family rather than a separate comparator.
**Comparator**: the deployed CVaR-VQE arm, `s24.d_harness.arm_vqe` at the S27 deployed settings
(alpha 0.18, T 0.5, layers 3, iters 80, lr 0.15, seeds 0 and 1), config DIS, whose per-target
realised m, `rmsd_vqe`, `rmsd_topm` and `rmsd_top75` are already in `s27/results/vqe_rows.jsonl`
(the artefact is reused, never recomputed; contract addendum 16).

## 2. The measurements
1. **The set statement.** Per target, the deployed arm's tail set against the top-m* prefix:
   exact set equality, |m(t) - m*|, and the Jaccard index. Reported as a distribution, not a mean.
2. **The point cloud.** Per target, |RMSD(control) - RMSD(deployed)|; the count within the
   built-chain input floor (S28-L18/R2: 0.006 A mean, 0.02 A on 12 of 126, 0.5 A worst case) and
   within 1e-6; the paired `ST.compare` with `ST.fmt` verbatim and the fold-clustered CI.
3. **The built chain, the reporting basis (rule 5).** The same contrast after
   `s12.instrument.project`, both sides projected IN THE SAME JOB from clouds built by the same
   code path, so the S28-L18 branch-flip floor cannot enter the contrast (S28-L43's discipline).
   Production's own chain row (`s27/results/chain_rows.jsonl :: DIS`, 3.2126) is quoted beside it
   and is NOT used as the comparator for the m* = 75 cell.
4. **Mechanism beside outcome (rule 18).** The realised m's distribution against m*, the
   correlation of m with chain length (native-free), and the m-ladder's own slope (how much RMSD
   one unit of m is worth), so "the quantum stage chooses m" is quantified in Angstroms per unit.

## 3. Falsifiers, registered before any contrast is computed
- **F-M6a (lane T's prediction, as relayed).** The deployed arm's emitted structure equals the
  fixed-profile control to within the built-chain input floor on **at least 120 of 126 targets**.
  Falsified if fewer than 120 agree at the floor for every m* in the registered sweep.
- **F-M6b (the endpoint claim, which is what the charter cares about).** The paired built-chain
  difference between the deployed arm and the best-agreeing fixed-profile control clears 0.7x its
  own MDE with the fold CI excluding zero. **If F-M6b does not fire, the deployed quantum stage
  contributes nothing at the endpoint that one target-independent number does not.**
- **F-M6c (the honest converse).** If the deployed arm is measurably BETTER than every fixed
  profile, the quantum stage contributes at the endpoint and this lane says so in the same entry,
  with the second seed as the replication.
**My registered prior.** F-M6a FAILS ON THE SET STATEMENT and PASSES ON THE STRUCTURE: exact set
equality will be rare (m varies with sd 6.8 around 74, and equals any single value on only a
handful of targets), but the emitted structures will nonetheless agree within the floor on most
targets, because the members entering or leaving the tail at rank ~70 to 80 carry ~1.3% of the
average's weight each. F-M6b does NOT fire (I expect |difference| well under 0.02 A on the chain,
inside its MDE). F-M6c does not fire. The interesting number will be the m-ladder slope: how
little RMSD the one surviving scalar is worth.

## 4. Multiplicity, cost, and what this does not claim
**Endpoint comparisons declared: 2** (the deployed arm against the registered-primary fixed
profile, on the point cloud and on the built chain), plus 4 sweep cells reported as a curve and
priced with `ST.best_of_k_within` if any is quoted as a best. The second seed replicates and is
not a separate comparison. Point cloud: zero compute (the artefact carries it). Built chain: 126
projections per arm in one governed job (`--agent S29D --tag CPU --est-ram 0.6`), per-target
checkpointed, one target probed first with its peak RSS quoted.
This experiment does NOT test whether a DIFFERENT quantum formulation could contribute (lane X's
configuration space is outside its scope), does not re-open the trainability result, and does not
claim anything about the quantum stage beyond the deployed candidate-index encoding at the
deployed (alpha, T). A null here is a statement about the production spine, which is exactly what
the contract's one hard constraint says S29 must build past.
