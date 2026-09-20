# S30 CONTRACT — the rules binding every lane

Inherited from S29's contract because each rule was paid for by an error. Where a rule exists
because of a specific incident, the incident is named — a rule whose cost is forgotten gets
dropped at the worst moment.

## Inherited non-negotiables

1. **Built-chain Cα RMSD is the endpoint.** Point-cloud numbers are intermediates and are labelled
   as such in the same sentence. *(S29's coordinator promoted a point-cloud ceiling into a headline
   and had to correct it — S29-L44.)*
2. **Below 0.7× MDE is not a result. 0.7–1.0× is not a demonstrated improvement.** This applies to
   our own positives as strictly as to anything else; it disqualified four of S29's.
3. **MDE is computed per comparison**, never per instrument. The project's old single 0.084 Å
   constant is wrong by up to 84× in both directions.
4. **Fold-clustered CIs govern**; i.i.d. CIs are secondary and reported beside them. Folds-same-sign
   is reported with every mean.
5. **Concentration is scored against a uniform-effect null**, never a raw drop-top threshold, which
   is not a valid test and has misfired here.
6. **Every controlled comparison uses a control in the operator's own space.** This is the project's
   most repeated error — three instances in two sprints, plus two more in S29 where a proposed
   control was the identity by algebra.
7. **Zero-information controls must be plausible, not uniform.** A worse measure is not an
   uninformative one.
8. **Per-target maxima over a family are order statistics**, priced with `best_of_k_within` and
   split-half transfer, never read as means.
9. **Never use the native to tune a deployable parameter.** ORACLE rows are labelled ORACLE in every
   sentence that uses them and kept in separate arms.
10. **Do not change benchmark60; do not regenerate the pinned folds or clusters.** Correcting an
    identity clustering once silently moved 13 benchmark targets and invalidated every fold model.
11. **The sealed benchmark is not spent.**
12. **Every experiment pre-registers its falsifier, committed before the first number exists.**

## Provenance and record-keeping

13. **An artefact path in prose is a claim, not a citation.** Before quoting the cost of resuming or
    reusing historical work, run `ls <path>` and `git log --all --oneline -- <path>`. *(Four paths
    named in this project's prose do not exist; one cost a mis-planned lane-week.)*
14. **Chronology is certified from git, not from anyone's word.** *(A coordinator claim that a
    formula predated the numbers it explained was false by 7.5 minutes.)*
15. **Corrections are annotated in place with the original wording left standing.** Nothing is
    silently edited away; a numbering gap is recorded, not renumbered.
16. **Commit by pathspec** — `git commit -m "..." -- path1 path2`, never `git add` first. Lanes
    share one git index and a failed commit otherwise leaves files staged for the next lane to
    sweep.
17. **Run `date` in the same command as a ledger or status append.** Guessing elapsed time has
    mis-stamped entries twice.

## Operational

18. **Check `s26/jobs/<name>.json` AND the live process list before launching.** A job name is not a
    lock: `jobrun.py` does not deduplicate, and S29 had three duplicate-job incidents in one night
    across three lanes with three different proximate causes.
19. **Temp files for concurrent writers carry the pid.** `os.replace` is atomic; a *shared*
    `<name>.tmp` is not, and an interleave inside it publishes atomically as corruption.
20. **Every expensive experiment is resumable and checkpointed.** Do not lose work to a resource
    kill.
21. **The box is shared with the user's own applications.** As of the S30 open they hold ~84% of
    RAM (Chrome 3.75 GB, VS Code 1.96 GB), leaving ~1.67 GB to the governor's 94% ceiling. Prefer
    many light jobs to few heavy ones, and re-measure headroom rather than assuming it.

## S30-specific

22. **The meter gates compute.** Every new cost function passes `s30/s30_D_meter.py` before
    receiving substantial VQE compute. A cost that moves in the wrong structural direction does not
    get a day of the box.
23. **A cost whose merits the meter cannot see is a finding about the meter.** Report it as such and
    extend the instrument; do not quietly route around it.
24. **State which S29 finding your hypothesis attacks and which it accepts as binding.** A
    hypothesis that engages neither the recognition failure nor the common-mode error carries a
    stated reason for bypassing both.
25. **Physical validity, energetic plausibility and native discrimination are three different
    things** and are never conflated. A cost that optimises the first two while ignoring the third
    reproduces exactly the failure mode S29 measured.
26. **Multiplicity is tracked sprint-wide.** With 8 lanes and many comparisons, a spurious 1× MDE
    positive is expected rather than surprising. Count the comparisons; hold positives to a bar that
    accounts for the count.
27. **The main team may not approve its own positive.** Anything that has not survived a genuine
    attempt to destroy it is labelled provisional in every document it appears in.
28. **Do not call something quantum merely because a quantum circuit produced it.** A quantum
    positive requires a matched classical control that fails to reproduce it, including a
    product-state/separable restriction to test whether entanglement is doing anything.
