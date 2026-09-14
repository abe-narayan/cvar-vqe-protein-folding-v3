# PREREG A1 -- DOES A qubit-ADAPT-GROWN ANSATZ CHANGE THE EMITTED STRUCTURE? (lane Q)

Filed 2026-09-13 09:00, before any endpoint number for any ADAPT arm exists. Not edited after;
addenda are appended. ENDPOINT experiment: the `--label` phase of `s26/q_adapt.py` reads the
native and is gated on "PHASE 0 SIGNED OFF" in `s26/LEDGER.md` (the script refuses otherwise).
The `--build` phase reads no native and no RMSD.

## 1. The question

Proposal A: replace the fixed 3-layer RY/CNOT-chain-plus-ring ansatz of the deployed CVaR-VQE
selector (`core/quantum.py` `StatevectorCircuit`, n = 7, P = 21) by a qubit-ADAPT-VQE ansatz
grown operator by operator on the same objective. Does the grown state change the emitted
structure, in which direction, and by how much against the comparison's own MDE?

## 2. What the record already knows (stated to the presenter whatever the outcome)

- The deployed Hamiltonian `E = _zrank(pool["sc"][top[:128]])` is the standardised rank ladder
  up to tie-averaging: worst deviation 4.06e-2 of an E range of 3.4371, 1.18% of range, 0 of 8
  targets exactly identical (`s25/results/q_gibbs.json`, `spectrum_target_independence`). Two
  trained states in the deployment, one per `(alpha, T)` cell of `VQE_LFO`: (1.0, 0.3) on folds
  0, 3, 4 (78 targets); (0.25, 0.3) on folds 1, 2 (48 targets) (`core/pipeline.py:113`).
- Selection is classical: the CVaR tail's support is a subset of an initial prefix of the
  energy order, 2,592 adversarial cells, 0 violations (`s25/results/q_verify.json`).
- The readout is insensitive to the state: the trained state sits 0.902 nats and 45% of mass
  from its Gibbs optimum (`s25/results/q_gibbs.json`), and the endpoint difference between the
  circuit and that optimum is -0.0302 A, SE 0.0450, 0.24x MDE, 43W/44L/39T, fold CI
  [-0.127, +0.053] (`s26/results/q_mde_reference.json`, contrast `s8:vqe_a1.0_T0.3-boltz_T0.3`).
- Deeper ansaetze and larger chi ordered nothing, five times over five sprints (S20 L1, L-B;
  S21 L16, L34; S22 L13 to L16; S23 L8).
- The optimiser trains: beats best-of-200 untrained draws at every T and closes 78 to 89% of
  the free-energy gap (`s25/results/q_gibbs.json`, `training_control`). No exponential
  plateau at n = 4..13 (`s25/results/q_plateau.json`). The DLA is being measured now (A2).
- New, from the pre-sign-off synthetic work (`s26/agentQ_FINDINGS.md` section 0): the
  deployed E is affine in the register index, index = rank, so exp(-E/T) factorises over the
  7 bits and the Gibbs state of the deployed Hamiltonian is a PRODUCT state: KL(Gibbs ||
  product of its marginals) = 7e-17 nats on the ideal ladder and 5.4e-5 nats on 1A13's real E
  (`s26/results/probe/1A13.json`, `product_diagnostics`). A single RY layer reaches it (KL to
  Gibbs 1e-4 after L-BFGS-B; `s26/results/probe/1A13.json` arm `adaptV_lbfgs_zrank_P7`), the
  deployed 21-parameter circuit does not (0.903 nats at 50 steps, 0.884 at 750). ADAPT at
  alpha = 1 selects no entangling operator and stops by its gradient criterion at P = 8 (V) or
  10 (L2).
- The prior expectation for A1 is therefore null. This experiment tests it rather than
  assumes it.

## 3. Hypothesis

H_A1: at matched parameter count (P = 21) a real-amplitude ansatz grown by qubit-ADAPT on the
deployed objective F = CVaR_alpha(E; p) - T H(p), with the deployed `(alpha, T)` per fold and
seed 0, emits a built chain whose CA-RMSD to the native differs from the fixed ansatz's by more
than the comparison's own MDE.

## 4. Arms (all per target, seed 0, the fold's `VQE_LFO` cell; every arm through the SAME readout)

    fixed_zrank_it50              core.quantum.run_cvar_vqe, the deployed selector. THE COMPARATOR.
    fixed_zrank_it750             the same, 750 Adam steps = the total Adam budget ADAPT spends
                                  (50 for the RY layer + 14 x 50). Budget control.
    adaptV_adam_best_zrank_P7/14/21    pool V (Tang et al. 2021, 2n-2 = 12 strings)
    adaptL2_adam_best_zrank_P7/14/21   pool L2 (all 1- and 2-local odd-Y strings, 91)
        re-optimiser: Adam lr 0.15, 50 steps per growth step, fresh state, BEST iterate
        (the deployed optimiser's arithmetic; the best-iterate rule is the one deviation,
        measured to be necessary in s26/agentQ_FINDINGS.md section 0.3); start = one RY layer,
        theta ~ N(0, 0.6^2) seed 0; growth stops at ||g_pool|| < 1e-3 or P = 21; ties in the
        operator selection broken by lowest pool index and counted.
    adaptV_lbfgs_zrank_P7/14/21, adaptL2_lbfgs_zrank_P7/14/21
        the same growth with L-BFGS-B re-optimisation (the ADAPT literature's choice;
        converges to |g| ~ 1e-7 where F is smooth). SECONDARY.
        If growth stops early, the final state stands in for the missing P and is flagged.
    gibbs_T                       exp(-E/T)/Z, the alpha = 1 analytic optimum, no circuit.
    uniform128                    uniform weights over the same 128: zero information in the
                                  operator's space.
    randH_fixed, randH_adaptL2    random weights softmax(c g), g ~ N(0,1) fixed per target,
                                  c bisected so the entropy equals fixed_zrank_it50's /
                                  adaptL2_adam_best_zrank_P21's: the matched random control.

PRIMARY contrasts (two): built chain, `adaptL2_adam_best_zrank_P21 - fixed_zrank_it50` and
`adaptV_adam_best_zrank_P21 - fixed_zrank_it50`.
SECONDARY: the same two on the selection basis; the lbfgs arms; P = 7 and P = 14; every
control against fixed_zrank_it50; the alpha = 1 (n = 78) and alpha = 0.25 (n = 48) subsets.

## 5. Falsifiers (both registered)

- "ADAPT helps" fires only if a PRIMARY contrast is negative, |effect| exceeds its own MDE,
  the fold-clustered CI excludes zero, 5 of 5 folds agree in sign, AND it replicates within
  its own CI on seed 1 and in reverse fold-processing order (`--seed 1`, `--order fold_rev`).
- "ADAPT is null" fires if both PRIMARY contrasts lie within +-0.5x their own MDE.
- Anything between is reported as underpowered or Type-M, never as a result.
- A positive on the alpha = 0.25 subset alone (n = 48, MDE about 1.6x larger) is reported as
  a subset finding with its own MDE and never promoted to the whole instrument.

## 6. Basis and readout

Built chain: `core.pipeline.average_weighted(Wo, block, p)` then `core.pipeline.project`
(ramah, lambda 0.3, multi-start, grad exact), scored as `core.pipeline.label` scores
`rmsd_q_synth`. The harness reproduces the production arm bit-for-bit: `ca` and `q_ca` equal
`bench_results/cache/464a0ddb5f283e04/<pdb>.json` at max |diff| 0.0 and the selection index is
equal on 1A13, 1A1P, 1CS9 and 2MK7 (`s26/results/probe/1A13.json`,
`s26/results/bitforbit/*.json`, key `cache_check`). Selection arm: `rr[o[consensus_medoid(
block, p)]]` = `rmsd_vqe_sel`, the s8-instrument readout, recorded beside it. The point cloud
`rmsd_q_avg` is recorded and is never compared with a structure.

## 7. Expected effect against the MDE (from stored artefacts; `s26/results/q_mde_reference.json`)

    built chain   rmsd_q_synth - rmsd_u_synth   -0.0135  SE 0.0342  MDE 0.0958  62W/64L
    built chain   rmsd_q_synth - rmsd_arm       +0.0133  SE 0.0181  MDE 0.0508
    selection     circuit - Boltzmann T=0.3     -0.0302  SE 0.0450  MDE 0.1261  0.24x
    selection     p-weighted - unweighted medoid -0.0308 SE 0.0590  MDE 0.1652

The fixed-vs-ADAPT paired SD is expected at or below the weighted-vs-uniform SD (both states
are trained toward the same objective), so the primary MDE is expected at or below 0.10 A on
the built chain; the run computes its own. Expected effect: on the 78 alpha = 1 targets ADAPT
converges to the Gibbs product state, so the contrast equals `gibbs_T - fixed`, whose
selection-basis analogue is -0.0302 A (0.24x MDE); on the 48 alpha = 0.25 targets the objective
is not the Gibbs functional and ADAPT-L2 reaches a lower F (-3.039 against -2.985 for the fixed
ansatz on the ideal ladder, `s26/agentQ_FINDINGS.md` 0.2), and the readout may move. Expected
|effect| < 0.05 A on both primaries: the null falsifier is expected to fire.

## 8. Memory and time

Peak RSS 0.338 GB for the full 16-arm harness on 1A13 (`s26/jobs_done/q_probe_1A13.json`),
0.379 GB fixed-only on 3 targets (`s26/jobs_done/q_probe_bitforbit3.json`). Wall 75.7 s for
1A13 under contention, 73.2 s of it projection (17 projections at about 4.3 s each,
`s26/results/probe/1A13.json`, `timings`). With the two randH arms, 19 projections, about 85 s
per target: 126 targets in about 3.0 h on one process, 1.6 h on two (two fit: 0.8 GB). Label
about 5 min (retrieve and score only). Stats 1 min. Replication on a positive: +3 h per rerun.
Agent-hours: 5 (null) to 11 (positive), mostly waiting. Launch:
`python s26/jobrun.py --agent Q --tag CPU --name a1_build --est-ram 0.5 -- python s26/q_adapt.py
--build --tag a1`; checkpoint per target under `s26/results/a1/<pdb>.json`, resumable.

## 9. Rule 0, six forks, each naming the alternative not taken

1. FUNCTIONAL. Per-target CA-RMSD of the emitted structure to the native, paired between
   arms. NOT TAKEN: the free energy, KL to Gibbs, entropy or objective value as the endpoint
   (all recorded as properties; an objective-level story is what was mistaken for a
   contribution in S8 to S24, S25 L5).
2. BASIS. The built chain through the production projection. NOT TAKEN: the point cloud
   (22.3% contracted, not a structure, S25 L4) and the AMBER-repaired emission (AMBER is not
   run; it would add a stage the contrast is not about).
3. READOUT. `average_weighted` then `project`, and `consensus_medoid(block, p)` for the
   selection arm, exactly as `core.pipeline.run_target` and `label`. NOT TAKEN: a
   probability-weighted CVaR-tail average (S23 L8 closed it in both regimes) or a uniform
   average over the realised tail (that is the classical top-m, S22 L16).
4. NORMALISATION. `E = _zrank` of the deployed score over the same 128: A1 asks about the
   ansatz, so the Hamiltonian is held at the deployed one. NOT TAKEN: raw z-score, asinh or
   soft compression of the score (A3's question; changing both at once would confound the
   ansatz with the Hamiltonian).
5. NULL. The fixed ansatz at its deployed 50 steps (comparator), fixed at 750 steps (matched
   Adam budget), `gibbs_T` (the analytic optimum, no circuit), `uniform128` (zero information
   in the operator's space), `randH_*` (matched entropy, no information about the candidates).
   NOT TAKEN: an initialisation mean or best-of-N untrained circuits as the ENDPOINT
   comparison (banned by `concentration-is-wrong-when-discrimination-binds`; best-of-N stays
   the objective-level training control S25 already ran).
6. LABEL. `s24.stats_lib.compare` paired at n = 126: SE, MDE = 2.8016 SE, iid CI beside the
   fold-clustered CI, W/L/T, median beside mean, drop-top-10 against the uniform-effect null,
   verdict. The medoid's `np.argmin` tie-break is kept as production has it (the bit-for-bit
   assertion needs it) and every selection-tie count is recorded; a secondary line averages the
   outcome over any tie set (`argmin_tied`). NOT TAKEN: marginal means, W/L alone, an iid CI
   alone, or a five-cluster CI excluding zero without the MDE (S24 D's own defect).

## 10. What I will not do

No cell of `VQE_LFO` is re-chosen. No ADAPT setting (eps, steps, pools, lr, start, optimiser)
is chosen by looking at an RMSD; all were fixed on the synthetic ideal ladder before this file.
The E variants belong to A3. The complex pool L2C is available in the code and is NOT an A1 arm.

## ADDENDUM (2026-09-13 21:50, after the run; the text above is unchanged)

Built 126/126 (`s26/results/a1/`, two shards killed by the host at 10:10 after 33 targets, one
governed process resumed from the checkpoints, `s26/jobs_done/a1_build.json`), labelled after
L33 (`s26/jobs_done/a1_label.json`), statistics in `s26/results/a1_stats.json`. The "ADAPT is
null" falsifier FIRED on both primaries (built chain: -0.0138 A at 0.23x MDE 0.0588, fold CI
[-0.071, +0.044]; -0.0222 A at 0.36x MDE 0.0608, fold CI [-0.085, +0.042]); "ADAPT helps" did
not fire; no replication owed. The expected effect (section 7: |effect| < 0.05 A, null
expected) held. Power: resolution 0.06 A on the built chain; underpowered below it. Ledger L68.

## ADDENDUM 2 (2026-09-13 22:35, BEFORE the run; robustness replication of a null, L77 scope)

A1 is null at the registered threshold, so no replication is owed by the contract; the user
extended the sprint (L77) and the record will ask. Run: `s26/q_adapt.py --build --tag a1s1
--seed 1 --order fold_rev --fixed-iters 50 --optimisers adam_best` (arms: fixed_zrank_it50,
adaptV/adaptL2 adam_best P7/14/21, gibbs_T, uniform128, randH_fixed, randH_adaptL2; 11
projections per target, about 50 s each, 126 targets, one governed process, checkpoints under
`s26/results/a1s1/`), then `--label --tag a1s1` and `--stats --tag a1s1`. Seed 1 changes the
fixed circuit's initial angles, the ADAPT RY layer's initial angles and the random-control
draws; reversed fold order changes only the processing order (each target is independent), so
the order half is procedural. What replicates: the two primaries, built chain. Prediction: both
inside +-0.5x their own MDE with fold CIs spanning zero, as on seed 0; the seed-0 point
estimates (-0.0138, -0.0222) lie inside the seed-1 iid CIs. Falsifier: either primary outside
+-0.5x MDE with a fold CI excluding zero on seed 1 (then the seed-0 null is not stable and the
discrepancy is reported as such). Launched after A3's build, as the coordinator ordered.

## ADDENDUM 3 (2026-09-14 04:12, after the addendum-2 run; ledger L139)

The addendum-2 falsifier FIRED: at seed 1 both primaries are outside +-0.5x MDE (0.71x, 0.79x)
with fold CIs excluding zero. The second prediction held (seed-0 points inside the seed-1 iid
CIs). A1 is NOT MEASURED on either seed; the seed difference sits in the fixed comparator
(3.228 -> 3.261 A), not in the ADAPT arms (3.214 -> 3.216). Verdict unchanged (L139).
