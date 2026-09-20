# S29 CONTRACT (read first; binds every S29 lane)

Sprint 29 of the peptide structure-prediction programme, repository
`C:\Users\abena\Protein-Folding-Algorithm`, branch `s26`, work under `s29/`. Python
`C:\Users\abena\miniforge_3\python.exe`. The authority is the user's charter, saved verbatim as
`s29/BRIEF.md`; this contract is its operational form plus the standing discipline inherited from
`s27/S28_CONTRACT.md` (all 12 non-negotiables and addenda 1 and 2 apply with `s29/` in place of
`s27/`; read that file too).

## The one hard constraint and the one endpoint
- The CVaR-VQE is the central scientific object: a parameterised quantum state, a Hamiltonian
  or observable set with a defensible meaning, a CVaR-family objective over a measured quantity,
  its optimisation causally connected to the final structure, and removing or randomising the
  quantum stage measurably degrades the result. A system whose quantum stage can be deleted
  without changing the output does not satisfy the charter (the production spine, where the
  CVaR tail is the classical top-m prefix to 1e-13, S28-L21, is exactly that; S29 must build past it).
- The endpoint is the mean built-chain CA RMSD on the 126-target instrument (production
  3.2126 A, `s27/results/chain_rows.jsonl :: DIS`; point cloud 3.0483). Point cloud is diagnostic.
  Every serious configuration reports its built-chain number, including the failures.

## Non-negotiables (unchanged from S28, restated)
1. Never open benchmark60. 2. Never regenerate or move a pinned artefact. 3. Production code is
frozen at `a15406c`; new work lives in `s29/` as `s29/s29_<lane>_<what>.py`, results
`s29/results/`, tests `tests/test_s29_<lane>.py`. 4. ORACLE is labelled in the same sentence,
every time; no native quantity chooses a deployable parameter. 5. Built chain is the reporting
basis. 6. `s24.stats_lib.compare`; fold-clustered CI decides; below 0.7x MDE is NOT A RESULT;
0.7 to 1.3x is the Type-M zone; grids priced with `best_of_k_within`; positives re-run on a
second seed. 7. Pre-register the falsifier before every run (`s29/PREREG_S29_<lane>.md`);
matched controls in the operator's own space; NaN-poison for every deployable operator.
8. Every python job over 200 MB or one minute goes through `python s26/jobrun.py --agent S29<X>
--tag CPU|AMBER|ESM|TEST --name <name> --est-ram <GB> -- python <script>`; one AMBER process at
a time; probe one target first and quote the peak RSS; checkpoint per target. The governor
(`s26/governor.py` v2.5) now runs the box at 92 to 94% RAM with a hard kill at 95.5%: the user
asked for 94 to 95% utilisation, so fill it, but keep every job resumable. 9. No slope is a
barren plateau or its absence; no quantum advantage claim; an exact simulator is never a cause.
10. Re-opening a closed question needs a genuinely new formulation, the reason the old closure
may not apply, a pre-registered falsifier and a meaningful control, stated up front in the
prereg with the sprint and ledger line that closed it. 11. Every number carries its artefact
path. 12. Ties never break by array order.

## S29-specific rules (from the charter)
13. Before any build, a lane states which of the charter's findings 1 to 11 its hypothesis
    attacks and which it accepts as binding; a hypothesis that ignores recognition (8) or the
    common-mode error (11) states why.
14. For every method imported from the literature: "what information does this method contain
    that the project does not?" If the answer is none, it is not imported. Rejections are
    recorded in `s29/lit/` with the reason.
15. For every serious quantum formulation, the nine questions of charter section 11 are answered
    in the ledger entry, and any accuracy positive is tested against the ten controls listed
    there (classical equivalent, diagonalised equivalent, permuted, random, matched budget,
    untrained circuit, simpler ansatz, order-statistic, product-state restriction, second seed).
16. Cheap and decisive before broad and expensive: a pre-registered probe on 12 targets before
    any 126-target run; a probe is never quoted as evidence for the instrument.
17. Multiplicity: every lane counts its comparisons in its ledger entries; the coordinator keeps
    the sprint total in `s29/STATE.md`; a positive at 1x MDE among many comparisons is held to a
    higher bar (state the max-over-K null).
18. Measurements of mechanism beside measurements of outcome: if a change helps, show the
    intermediate quantity that moved.
19. The cost-RMSD meter (`s29/s29_D_cost_audit.py`, lane D) reports four numbers for any cost
    function: the ladder Spearman, the gradient cosine at production, the native percentile in
    the pool, the preference for the ORACLE structure with the pool-member control. Every new
    objective goes through it before any endpoint run; the meter is ORACLE and never tunes.

## Ledger, state, status
- `s29/LEDGER.md` (append-only; `## S29-L<n> -- TITLE (date time, lane)`; number from the tail
  in ONE python process; a colliding number takes a b-suffix). Every entry: question, falsifier,
  result with `ST.fmt` verbatim where a contrast is claimed, verdict, artefact path.
- `s29/STATE.md` (coordinator): leading hypothesis, strongest objection, what would falsify it,
  each agent's assignment, what is closed, the comparison count. Updated as results arrive.
- `s29/STATUS.md`: two lines per lane per hour. `s29/s29_<lane>_FINDINGS.md`: the lane's
  findings in the standing format. `s29/lit/<lane>_<topic>.md`: literature notes (paper, the
  equation that matters, the assumption checked against this pool, kept or rejected and why).
- Commit your own files early and often on branch `s26` with the single trailer line
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`. Never `git add -A`.
- Every lane finishes each turn with two lines: what stands, what is next (or what fell).
