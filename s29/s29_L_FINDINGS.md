# S29 LANE L -- FINDINGS (literature, permanent role)

Lane L read for the whole sprint and ran no experiment. Eight topics, ~60 primary papers read for
their constructions rather than their abstracts, in `s29/lit/L_1` to `L_8` with
`s29/lit/L_INDEX.md` as the one-line-per-paper index. Ledger: S29-L1, L8, L12, L13, L14, L16,
L19, L31, and the correction L42.

Format follows the S12-to-S25 convention: DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS /
REFUTED / OPEN, then what damaged my expectations, then what I did not do and why.

---

## 1. DEMONSTRATED (from the literature, checked against this instrument)

**L1. No published model-quality method has ever been trained or benchmarked below about 40 to 50
residues.** ProQ3 filters out every target under 50 residues by name; VoroMQA's learning set is
chains longer than 99; DeepAccNet trains on 50 to 300; QMEANDisCo quotes ~0.12 expected error at
40 residues as its small-model limit. The four signal classes they reduce to are each absent at
9 to 16 aa (packing/burial), already this project's own axis (sequence-to-local-structure
agreement), falsified here by measurement (consensus, against a 68% common-mode pool), or measured
worse than random here (single-point physics). Charter finding 8 is the field's position at this
length, not a defect of the S27 library. (S29-L1)

**L2. AlphaFold2's pLDDT has no within-target ranking skill on 588 peptides of 10 to 40 aa**
(McDonald et al., Structure 31:111, 2023): "the lowest RMSD structures failed to correlate with
lowest pLDDT ranked structures"; taking rank-1 instead of best-of-5 costs 0.2 to 1.1 A. The best
confidence signal in the field, from the model that solved the protein case, fails in-band at
peptide length. External confirmation of finding 8 on a different instrument. (S29-L1)

**L3. S23 L9's error identity IS the Krogh-Vedelsby ambiguity decomposition**, and the
Ueda-Nakano form supplies the (1 - 1/M) coefficient on the covariance that does not decay. Finding
11 is a universal law, not a property of this pool. Derived consequence (arithmetic on
`s23/results/errdecomp.json`, caveated in the note): an infinite pool of the same kind returns
3.040 A against the shipped 3.0483, so "average more members" is worth <= ~0.008 A. (S29-L8)

**L4. Every literature method that genuinely breaks correlated error needs one of four inputs and
this instrument has none:** trainable members (negative correlation learning), a control variate
whose mean is KNOWN (control variates -- the known mean here is the native), samples of the
high-fidelity quantity (multifidelity Monte Carlo, which is unbiased precisely because it keeps an
anchor), or a known bias ratio along a shared direction (Richardson-style two-source
extrapolation, which is H1's family). (S29-L8)

**L5. The perception-distortion theorem makes finding 8 a necessity.** Blau & Michaeli (CVPR 2018)
Theorem 3: for ANY distortion measure, the distortion-optimal estimator's output DISTRIBUTION must
diverge from the distribution of real signals, most steeply at the low-distortion end where this
project operates. Every scorer in the S27 library is a realism measure, so all of them must
disprefer the RMSD-optimal answer. The 25.8% averaging contraction is Jensen's inequality, not a
distogram bias; post-hoc calibration is closed by argmax/median invariance, which explains S25 L2
rather than merely repeating it. (S29-L12)

**L6. The project's set-equality theorem is eq (12) of Barkoutsos et al. (2020)** -- CVaR is
defined on SORTED samples and the paper scopes itself to "classical optimization problems, which
yield diagonal Hamiltonians". Cite it; it is a definition, not a discovery. The same paper gives a
falsifier nobody here had used: CVaR's global-minimiser set is {theta : overlap >= alpha}, large
and flat, so any accuracy change attributed to CVaR optimisation must be shown not to be a
tie-break inside that set. (S29-L13)

**L7. Non-classicality requires BOTH non-commuting terms AND a target that is not an eigenvector.**
Gibbs/thermal states and quantum Boltzmann machines satisfy the second; S28's non-diagonal
Hamiltonian satisfied the first with a degenerate rank-one off-diagonal. The project has never had
both. And arXiv:2312.09121 (Nat Commun 16:7907, 2025) makes "our circuit trains well" evidence
FOR classical simulability, so the charter's ten controls are the central test, not a formality.
(S29-L13)

**L8. There is no published ceiling for native-free peptide prediction at 9 to 16 residues.** The
published selected-answer numbers at 9 to 25 aa are PEP-FOLD 2.6 A (25 NMR peptides) and APPTEST
1.96 A (42 peptides), neither like-for-like: different benchmark composition (56% of FAIL18 is
amyloid/lasso classes those sets exclude), best-of-N versus single-answer reporting, and every
published method selects within an ensemble its own energy generated. The comparison splits --
our GENERATION (ORACLE 1.71 to 2.31 A) is inside the published band; SELECTION is the whole 0.9 to
1.5 A gap. (S29-L14)

**L9. Nobody has conditioned on a native-free realism statistic before measuring accuracy
ordering.** The closest methodological paper (Hamelryck et al., PMC2677743) documents the
confounds -- for 139 of 149 decoy sets the native is trivially discriminable from individual
energy terms -- and explicitly declines the matching fix. Lane D's measurement is novel. The
design itself is published and human-validated (PIRM 2018), and its own result is the warning:
the realism index tracks human judgement at Spearman 0.83 BETWEEN bands and unreliably WITHIN
them. (S29-L16)

**L10. The in-band signal that survives a band has a closed form: the partial correlation.** For
jointly Gaussian (S, Y, R) the conditional correlation given R is constant across the band and
equals rho_SY.R = (rho_SY - rho_SR rho_RY)/sqrt((1 - rho_SR^2)(1 - rho_RY^2)), exactly zero iff a
scorer's association with accuracy is fully mediated by realism. Band width interpolates between
the global correlation and the partial correlation, so report skill as a curve over widths.
(S29-L16)

**L11. THE PER-TARGET SIGN IS AN INCIDENTAL PARAMETER (Neyman & Scott 1948), so it is not
estimable from other targets' answers as a matter of statistical theory, not of model capacity.**
Pooled ML is inconsistent for it; the standard remedy (conditional/fixed-effects likelihood =
the within-group pairwise design) ELIMINATES it rather than estimating it. S14's 0.986-within /
0.600-across with a -0.0005 overfitting gap and a learning curve that declines with more targets
is that theorem measured. The only two escapes the literature names: replication within the
instance, or a covariate observed at inference. **This is the sprint's unified finding with a
proof behind it, and I regard it as lane L's most durable output.** (S29-L31)

**L12. The field's per-target conditioner is MSA depth, and it is structurally absent at this
length.** AF2 degrades below MSA depth ~30; a 13-mer has no family, and any hits are the
fragment's parent proteins, i.e. leakage. Recycling is OPTIMISATION, not information -- it feeds
back only the network's own output. This is the cleanest statement of why the peptide case is
harder than the protein case. (S29-L31)

## 2. ORACLE DIAGNOSTIC / arithmetic (labelled, not measurement)

**L13. Pricing a per-target sign predictor**, using lane T's formula
gain = RMSD(1 - sqrt(1 - rho^2 (2q-1)^2)) at rho = 0.37 and RMSD_prod 3.2126, verified to
reproduce T's own 2.985 A at q = 1: q = 0.60 -> 0.009 A, 0.70 -> 0.035, 0.80 -> 0.080, 0.90 ->
0.144, 1.00 -> 0.228. The gain is QUADRATIC in (2q-1). **I published this table four times too
large in a first draft and corrected it in the same entry** (S29-L31). The correction changes the
conclusion: a genuinely good sign classifier buys less than this instrument's noise, and the
ceiling of the whole direction is 2.985 A.

## 3. HYPOTHESIS (mine, unmeasured, flagged as such)

**L14. The entropy term may be compactness-like**, because what conformational entropy tracks at
peptide length is basin width, and compact conformations have narrower basins. If so the
free-energy route is not orthogonal to the realism axis, and both of the sprint's live ideas act
on the same axis in opposite directions. Lane T is measuring the owned-channel version of this
(32 channels x 500 members x 126 pools, native-free Spearman with member Rg); its four-target
smoke test (LEG_compactness +0.956, RG_LAW +0.904, LEG_solvation +0.766) is consistent with it.
(S29-L19, S29-L42)

**L15. A compactness-loaded realism band would produce a self-fulfilling null.** S14's flip
diagnostic found per-target in-band skill correlating +0.909 with the native's z-scored Rg, so the
in-band axis IS compactness; every realism statistic in the library is compactness-like. Banding
on one would crush rho_SY.R to zero BY CONSTRUCTION. The fix is to measure rho(R, Rg) for the
chosen band statistic and report it beside the result. (S29-L19)

## 4. REFUTED / corrected (my own)

**L16. "The S8 free-energy stage is committed and resumable" -- FALSE, my error.** All seven S8
paths are absent from disk and from all of git history; `s8.relax` is not importable. What
survives is prose in `docs/FINDINGS.md` section B. The item is a REBUILD from a prose spec, not a
resume, and the corrected gate is that the stage must be BUILT before any diagnostic can run.
(S29-L42, after lane T's S29-L41; confirmed independently by me.)

**L17. I nearly re-imported a closed route.** I drafted PEP-FOLD's retrieval key (retrieve by
predicted local conformation) as the strongest unexploited peptide-literature lead. It is closed
here twice -- S13's torsion-bin key (0.517 on FAIL18 against a 0.562 majority baseline) and the
22-key / 24-arm screen (no key beats BLOSUM on the chain; structural keys worsen pool-best
1.711 -> 2.161 and the ORACLE ceiling 1.994 -> 2.542). Caught before committing. (S29-L14)

## 5. OPEN

- Whether an independent second generator is obtainable. Everything in the quasi-single-model QA
  family (CASP15's EMA winners) and escape route (1) of L11 depends on it, and its value is the
  independence, not the scoring rule. S24 L2/L3 is the measured blocker, not a proof of absence.
- Whether a superposition-free local objective (an lDDT-like quantity) behaves differently in the
  cost meter. The CASP assessors say plainly that superposition-free local measures are the
  learnable ones and GDT-TS is the hard one; our objective and endpoint are both
  superposition-dependent. Cheap for lane D and never asked here.
- The anisotropic reading of pool dispersion. Spread is a MAGNITUDE statistic and is
  reflection-symmetric, so it cannot reach the sign; only the orientation of its principal axes
  could. Not what the existing memory note measured.

---

## 6. WHAT S30 SHOULD READ -- and the lesson that belongs beside the list

### The reading list, shortest useful form

If S30 reads five things, read these, in this order.

1. **Neyman J, Scott EL, Econometrica 16:1-32 (1948).** Why the per-target sign cannot be learned
   from other targets. It converts this project's central empirical closure into a structural one
   and tells you the only two escapes. `s29/lit/L_8_conditioning.md` section (b).
2. **Blau Y, Michaeli T, CVPR 2018 (arXiv:1711.06077), Theorems 1 and 3.** Why a realism scorer
   must disprefer the RMSD-optimal answer, for any distortion measure, most steeply where we
   operate. `s29/lit/L_3_decision_theory.md`.
3. **Brown, Wyatt & Tino, JMLR 6:1621-1650 (2005), eqs (9) and (10).** The ambiguity and
   bias-variance-covariance decompositions -- i.e. this project's own S23 L9 identity, and the
   (1 - 1/M) coefficient that caps every aggregation route. `s29/lit/L_2_correlated_error.md`.
4. **McDonald EF et al., Structure 31:111-119 (2023).** The field's only ceiling-like statement at
   peptide length, and it is a ceiling on SELECTION. `s29/lit/L_1_native_free_qa.md` section 4.1.
5. **Cerezo M et al., Nat Commun 16:7907 (2025), arXiv:2312.09121.** Why trainability is now
   evidence FOR classical simulability, and therefore how any quantum claim must be positioned.
   `s29/lit/L_4_quantum.md` section 2.4.

Everything else, with the information test applied per paper, is in `s29/lit/L_INDEX.md`.

### The lesson that belongs beside it

**A FINDINGS paragraph asserting that a module is committed is not evidence that it is.**

This is the second instance of the failure mode in this project. The first is S26's presentation
file, recorded in project memory as "nowhere on disk or in git". The second is mine: I read
`docs/FINDINGS.md` section B -- which names a module, its tests, four artefacts and six columns,
in the present tense -- and reported the stage as "committed and resumable", which propagated into
a coordinator assignment. Neither the code nor its one completed target has ever existed in git
history. Lane T caught it by running an existence check as its first step, and said it did so only
because last sprint's memory told it to.

The specific trap is that detailed, internally consistent prose about a stage that was real in a
live session reads exactly like verification. The generalisation for S30:

- **An artefact path in prose is a CLAIM, not a citation.** Contract rule 11 ("every number
  carries its artefact path") is worth precisely what an existence check on that path is worth.
- **Run the check during the reading, not after it.** `git log --all --oneline -- <path>` costs
  one second and is the only thing that distinguishes a committed module from a remembered one.
- **The risk concentrates exactly where it hurts most**: on the one item that survived a
  closure, because that is the item the next sprint will act on. Both instances of this failure
  mode were on a surviving actionable item, not on a closed one.

## 7. What damaged my own expectations

- I expected to find at least one importable native-free scorer in the QA literature. There is
  none for this length; the entire field starts at 40 to 50 residues. The negative was worth more
  than the search would have been.
- I expected the ensemble literature to offer a way to break common-mode error. It offers four,
  and this instrument has the inputs for none of them.
- I expected the peptide literature to hold an unexploited architectural lead. Its central claim
  was already tested and closed here, twice.
- I expected a sign predictor to be worth more than it is. The quadratic dependence on (2q-1)
  makes even a good classifier worth less than the noise, and I had the table wrong by 4x before
  I checked it.

## 8. What I did not do, and why

- No experiments, by brief. Every number in my files is cited to a paper or to a prior sprint's
  artefact; the only arithmetic I produced (the <= 0.008 A averaging bound and the sign price
  table) is labelled as arithmetic on existing artefacts, with its caveats, and the sign table was
  cross-checked against lane T's independently published 2.985 A.
- I did not run the rho(R, Rg) check, the DLA measurement, the free-energy orthogonality
  diagnostic or the sign regression, though I recommended all four. They belong to lanes with
  experimental remits, and three of the four were picked up (T took the DLA and the orthogonality
  question; D took the band design).
- I did not verify every artefact path I cited from prior sprints' prose. I verified the one that
  mattered only after lane T found it missing, which is the wrong order and is the subject of
  section 6.
