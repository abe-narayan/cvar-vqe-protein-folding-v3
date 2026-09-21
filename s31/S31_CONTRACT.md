# S31 CONTRACT

Every rule is annotated with the incident that paid for it. Rules without an incident are charter
text. **A rule you cannot cite an incident for is a preference, not a contract.**

---

## I. Inherited non-negotiables

**1. The endpoint is the mean built-chain Cα RMSD on the 126 `tuning126` targets — 3.2105 Å.**
The CA point cloud (3.0483) and the *set mean* (3.5507) are different objects. State the basis in
the same sentence as the number.
*Paid for by:* S29 quoted the architectural ceiling on the cloud basis when the endpoint is the
chain; S30 shipped three live production values before lane V forced the reconciliation.

**2. The cloud→chain transfer for a correction is 0.92, measured on two arms.**
Do **not** use the 1.16 from `operator-consumes-set-mean` — that maps *set mean* to output, and its
own inputs give 4.21 Å against production's actual 3.21.
*Paid for by:* my S30 transfer of a cloud delta with the wrong law; caught by lane P.

**3. MDE = 2.8016 × SE, per comparison. Below 0.7× is not a result; 0.7–1.0× is NOT MEASURED.**
Fold-clustered CIs on the pinned folds. `stats_lib._verdict` refuses a verdict without folds.
*Paid for by:* S29's sprint-wide 0.084 Å constant, wrong by up to 84× in both directions.

**4. A number below its own MDE may support a NULL. It may never support a presence.**
*Paid for by:* I twice used a sub-MDE number as evidence of absence in S30 — once at 0.22× MDE, in
the same document that states the 0.7× rule two pages earlier.

**5. `benchmark60` is sealed; its single pre-registered pass is spent.** Folds and clusters are
pinned and are never regenerated. The native never tunes a deployable parameter.
*Paid for by:* correcting the identity clustering once silently moved 13 benchmark targets and
invalidated every fold model.

**6. ORACLE labels travel inside the sentence carrying the number.**
*Paid for by:* S30 published the `coh` gate as "native-free" when both its arguments require the
native; two of six appearances of its price dropped the NOT-DEPLOYABLE label, and in exactly the two
places framed as *what you get if you pass the gate*.

**7. A control must match the operator's own space.** The project's most repeated error.
*Paid for by:* S30's "retrieval is exonerated" compared an ORACLE best-of-pool CA-cloud measurement
against a record about the emitted built chain, at power 0.05 — while §6.1 of the same report
claimed there were no instances that sprint.

**8. A statistic must be matched to its null.** Five instances across S29–S30: a right-skewed
statistic's median read against a mean's baseline; `abs(signed mean)` scored against a null for
`mean |cos|`; a 0/0.5/1 indicator whose median is uninformative by construction; `stats_lib.compare`
being lower-is-better and inverting the verdict string for preference rates; and lane P hitting that
same inversion *in a lane that had just read the entry about it*.
**Quote the gate, never `.verdict`.**

**9. A zero-information control must be plausible, not uniform.** For torsions, a constant α-helix;
uniform-on-the-torus is a *worse* measure, not an uninformative one.

**10. A control that is itself a random draw needs its own draw distribution, and the number of
draws scales inversely with the effect size.**
*Paid for by:* S29 published the **maximum of its own eight** possible draws; eight separated the
two cases, one did not. Relative draw noise was 59% on the chain against 27% on the cloud.

**11. A per-target maximum over K variants is mostly best-of-K.** Any sweep needs a split-half
transfer arm; W/L cannot diagnose it.

**12. A stratum defined by the outcome cannot measure the thing that defines it.**
*Paid for by:* `FAIL18` is defined by the filter's own recall (`s12/instrument.py:271-278`), and
S30-L2's headline was 49% forced arithmetic. **A matched control in the right space does not rescue
a stratum defined by the outcome.**

---

## II. Provenance

**13. Historical artefacts are immutable.** Corrections are annotated in place with the original
wording left standing. Do not delete inconvenient results or overwrite prior ledgers.

**14. Prose is not evidence of code.** Six instances, the sixth in *shipped code*
(`core/pipeline.py:821` still asserts a number S25 withdrew). Run `ls` and read the file.

**15. Check the job, not just the file.** A file being written by a live job is not evidence about
that job; missing, unfinished and crashed look identical to `ls`.
*Paid for by:* S30's synthesis lane declared a **completed** 126-row result "unresolved at 112/126"
because it read the rows mid-write, and nearly invalidated the sprint's headline.

**16. `x.get(k) or 0` manufactures data.** It turns a missing key into a confident 0.0.
**Repeated identical values across different questions is the signature of a null input** — three
the same is suspicious, sixteen the same is a bug.
*Paid for by:* lane R's 0.500/0.500/0.500 cell, and my own 32-cost sweep printing sixteen fabricated
`native_pctile = 0.0000` three hours after recording lane R's case.

**17. Every meaningful number carries:** artefact path, experiment ID, config, seed, endpoint basis,
native-free-or-ORACLE, and its comparison baseline.

**18. Ledger numbers collide.** Three rounds in one S30 afternoon. Take the next free number,
re-check immediately before writing, renumber the later-stamped entry on a collision, annotate in
place.

---

## III. Method

**19. Pre-register the falsifier in a committed file before the number exists.**
*Earned:* four of S30's ten registered falsifiers fired against the lane that wrote them. **A
pre-registration that only ever confirms is decoration.**

**20. Work the ladder:** derivation → toy → small real subset → fold-held-out → 126 → deployment.
Do not jump from a toy to the endpoint.

**21. An implied-endpoint conversion from ρ is an upper bound attained only by a perfectly
calibrated correction. Always also APPLY the correction and measure.**
*Paid for by:* S30's implied 3.0338 against an applied 3.0519 — the conversion flattered a no-op.

**22. A pre-check that can only remove your own excuse is worth running before the contrast, not
after.**
*Earned:* lane G's chirality occupancy pre-check refuted its own lane's registered mechanism and
**strengthened** the conclusion by doing so — it converted a quiet null into a real negative.

**23. Multiplicity is tracked sprint-wide, written to as comparisons are emitted.** S30's total was
only ever a lower bound in the hundreds because each lane counted its own.

**24. The main team may not approve its own positive.** Two independent adversaries found **41
defects** in the S30 report, six severe, four of which inverted a headline. Budget for this.

**25. A withdrawal is a claim too, and can be wrong in the same ways.**
*Paid for by:* two lanes independently talked themselves out of the divergence claim — one on a null
about a different *outcome*, one on an argument about a different *mechanism* — and it was then
confirmed at 3.56× MDE. **Two revisions agreeing is not evidence when each rests on a different
object than the one under test.**

---

## IV. S31-specific

**26. "Make H non-diagonal" is not a result.** S30 closed it: a graph Hamiltonian gave a nearly
rank-one hopping term with negligible gradient contribution. A non-diagonal `H` must be physically
meaningful and must answer the charter's §10 twelve questions.

**27. Determine classical reducibility explicitly (charter §11).** For a diagonal `H`, CVaR reduces
to classical tail selection. **A classically reducible construction may still be useful — say so
plainly — but must not be described as a quantum mechanism.** Give the classical algorithm and its
cost.

**28. State which source a claimed new observable draws on.** S30's enumeration collapses to
**sequence** and **library**, with physics an *operator* on either. A new channel must name its
target argument.

**29. Correlating with the error is not the test. Incoherence with the pool's common mode is.**
At identical R² a coherent corrector emits +0.0554 worse and an i.i.d. one −0.2466 better.

**30. Do not claim an information-theoretic impossibility without proving one.** S30's honest form
is *"the extractable part is the wrong part, by a theorem about which part is identifiable"* —
which is stronger than an exhausted search and weaker than a mutual-information zero. Keep theorem,
empirical estimate and ORACLE quantity in separate categories.

**31. No sub-0.01 Å chain claim is valid until the projection seed is pinned.** The spread across
instruments is 0.0107 Å against this project's one confirmed effect of 0.0221 Å.

**32. Do not engineer a fragile 2.98.** A ceiling argument that redirects the next three sprints is
worth more, and the charter says so twice.
