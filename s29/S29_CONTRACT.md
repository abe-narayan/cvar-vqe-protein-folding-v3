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

## Addendum 1 (2026-09-20 00:05, coordinator; binds every lane from now)
20. THE SHRINK RULE (from lane T's theorem 2, S29-L7). The meter's gradient cosine is GAMEABLE:
    shrinking an objective's target map toward typicality, m -> tau + s(m - tau), drives the
    cosine positive with ZERO information added while moving the emitted structure toward the
    typical map. Therefore: every candidate objective that reports a gain in meter number 2
    (the cosine) MUST report, in the same entry, (a) its implied shrink s (the scale of its
    target map's deviation from the length-matched typical map, measured against the shipped
    map's), (b) meter number 3 (the native's percentile in the pool), which the shrink moves
    the wrong way, and (c) the emitted structure's mean virtual bond and Rg. A cosine gain
    without those three is not readable and lane D will veto it on sight.
21. THE MARGINAL CLASS IS BOUNDED (same theorem). No objective built from the per-pair
    marginals of the shipped posterior -- separable or not, any functional form, any
    normalisation -- has an expected gradient cosine whose magnitude is set by the native's
    deviation from typical; the expectation is second order in the pool's own spread and is
    exactly zero when the posterior's median map and production deviate from typical alike.
    A lane proposing a new objective INSIDE that class states what it expects to buy other
    than local informativeness (ordering, variance, non-contraction), and does not claim the
    gradient. A lane proposing to leave the class names the channel that sees the native's
    deviation from typical and prices its error coherence (S19).

## Addendum 2 (2026-09-20 00:22, coordinator; from lane L's topic 3 and topic 4)
22. THE PERCEPTION-DISTORTION BOUND (Blau & Michaeli 2018, theorem 3; lane L, S29-L12). For ANY
    distortion measure, the distortion-optimal estimator's output DISTRIBUTION must diverge from
    the distribution of real signals, most steeply at the low-distortion end. Every scorer in the
    S27 library is a realism measure, so every one of them MUST disprefer the RMSD-optimal
    answer: S28-L48 is that theorem measured, and the scorers are not defective. CONSEQUENCE: no
    lane may propose "find a native-free scorer that prefers the near-native structure" as a
    route; it is forbidden in advance. The one open version is ordering WITHIN a matched-realism
    band, which lane D is measuring. Any future recognition claim states which side of this
    bound it sits on.
23. THE AVERAGING BOUND (Krogh-Vedelsby and Ueda-Nakano; lane L, S29-L8). S23 L9's error
    identity IS the ambiguity decomposition, and the (1 - 1/M) coefficient gives an infinite
    pool of the same kind at 3.040 A against the shipped 3.0483: enlarging or re-weighting the
    retained set is worth at most about 0.008 A. No lane spends a 126-target run on "average
    more or differently" without stating why it escapes this bound.
24. THE FLAT-MINIMISER CLAUSE (Barkoutsos et al.; lane L, S29-L13). CVaR_alpha's global-minimiser
    set is {theta : overlap with the best state >= alpha} and is large and flat. Any accuracy
    change attributed to CVaR optimisation must be shown NOT to be a tie-break inside that set.
25. CONTRACTION IS JENSEN, NOT THE POSTERIOR (lane L, S29-L12; lane D, S29-L10). The Bayes
    estimator of a symmetrically over-confident posterior is not contracted; the 22 to 26%
    contraction is Jensen's inequality on the coordinate average, and the shipped cost's
    gradient EXPANDS. No lane may claim a contraction mechanism at the posterior or gradient
    level; the arithmetic of averaging is the mechanism.

## Addendum 3 (2026-09-20 00:52, coordinator; from lane P's registered clause)
26. THE BOUND DOES NOT LICENSE DISMISSING A MEASUREMENT. Lane T's achievable bound (S29-L23)
    now predicts a null for most remaining arms. That prediction is a hypothesis with a
    load-bearing assumption (B2), not a verdict. Therefore, for every arm already
    pre-registered before the bound existed: the falsifier is honoured exactly as written, no
    prior is revised, no threshold is moved, and no arm is dropped because the bound says it is
    null. IF AN ARM CLEARS ITS FALSIFIER IT IS REPORTED AS EVIDENCE AGAINST THE BOUND, with the
    second seed and the tie-key replication run FIRST, rather than explained away. A bound that
    forbids a measured effect is the bound's problem. Lane P registered this clause itself
    before its numbers existed and it is adopted sprint-wide.
27. CHRONOLOGY IS PART OF THE EVIDENCE. Any entry claiming "this was predicted in advance" must
    carry the timestamps of the prediction and of the measurement, so a reader can check that
    the prior was not moved after the fact. Lane P's addendum 5 is the model: falsifiers and
    null prior registered 00:10, its own mechanism prediction 00:37, lane T's bound 00:43.

## Addendum 4 (2026-09-20 01:28, coordinator; from lane D's S29-L34)
28. COMMIT BY PATHSPEC, NEVER BY INDEX. Eight lanes share one working tree and therefore ONE git
    index. `git add` followed by `git commit` is not atomic: if the commit fails (a transient
    object-lock from a concurrent lane is common here) the files stay STAGED, and the next lane
    to commit sweeps them into ITS commit under ITS message. This has now happened at least
    three times tonight (S29-L10 committed as "S29-L8", S29-L28 as "S29-L27", and S29-L33's
    band entry and its three artefacts inside lane O's commit e2e49109).
    THE RULE: commit with an explicit pathspec in ONE command --
        git commit -m "<message>" -- s29/path/one s29/path/two
    which commits exactly those paths and does not touch the shared index. Do not `git add`
    first. If a commit fails, re-run the same one command; nothing is left staged for another
    lane to inherit. The contract's existing "never git add -A" does not cover this case,
    because the hazard is another lane committing what YOU staged.
29. THE LEDGER HEADING IS AUTHORITATIVE, NOT THE COMMIT MESSAGE. Where the two disagree (the
    three cases above), the entry's own heading and artefact list decide, and the discrepancy is
    recorded in a provenance note rather than corrected by rewriting history. The final report
    traces numbers to LEDGER ENTRIES and artefact paths, with commits as corroboration.
