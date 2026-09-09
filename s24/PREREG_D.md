# PREREG — WORKSTREAM D, SPRINT 24

Written before the run in every case. Rule 0: six fork axes, each NAMING THE ALTERNATIVE NOT
TAKEN. Sent to the coordinator before execution.

---

## D0 — SET-EQUALITY DETERMINATION (`d_setequality_proof.py`) — NOT DIRECTIONAL

Registered as a **property test**, not an experiment: it reads no RMSD, touches no candidate
pool, and has no outcome that can favour a hypothesis. Rule 0's fork enumeration does not
apply and is not claimed. What is registered instead is a **falsifier**:

> Exhibit one `(E, p, alpha)` triple whose realised CVaR tail support is NOT a subset of an
> initial prefix of the energy order — i.e. a state carrying tail mass whose energy exceeds
> that of some state carrying none.

**REV1 (superseded).** Asserted the stronger claim — support is always an *initial prefix* —
and reported 0/2916. **That assertion was dead**: every probability family carried an epsilon
floor (`+1e-30`, `1e-18`, `+1e-6`), so `p > 0` held everywhere and the assertion could not
fire. Corrected by the coordinator; the over-claim is retained as a *measured statistic*.

**REV2 (standing).** Asserts subset-hood plus "every hole is exactly a zero-probability
state". Falsifier did not fire: 0/3888 on both. Prefix-hood violated on 1164/3888 (29.9%),
which is the coordinator's result reproduced. Full-support cells: 2268/2268 exact equality.

**Standing statement (coordinator's wording, used verbatim everywhere):** *The realised CVaR
tail's support is always a SUBSET of an initial prefix of the energy order, and equals that
prefix exactly when every state in the prefix carries positive probability. p_theta can delete
a member; it can never add one outside the classical top-m.*

---

## D1 — HAMILTONIAN-DISAGREEMENT SAMPLING (`d_hamiltonians.py`) — DIRECTIONAL

Directive §28, promoted by the coordinator to Lane D's highest-value item after L3 showed that
selection by the shipped score makes candidate sets parallel regardless of provenance
(+0.9432 vs a +0.9330 within-source control). If provenance is worth nothing and the selection
rule is fixed by D0's theorem, then **the functional is the only remaining lever**, and Legacy
and AMBER are the only two genuinely different functionals available.

### The hypothesis, stated so it can fail

**H_D1.** Partitioning candidates by *Legacy/AMBER disagreement* isolates subsets whose emitted
coordinate average has a **bias direction materially non-parallel to the incumbent's**, i.e. a
bias cosine against the incumbent below the ~0.65 bar L2 set, measured against the +0.9330
within-source control.

**Pessimistic prior, stated in advance.** I expect this to FAIL the bar. Both energies are
full-register in torsion space, both are dominated by compactness and sterics on real protein
windows, and L2/L3 have already shown that anything selected to be *good* ends up parallel.
The most likely outcome is cosines in the 0.85–0.95 band — i.e. disagreement partitions differ
in QUALITY but not in DIRECTION. If that is what comes back, it is the result, and it closes
the functional lever the same way L3 closed the provenance lever.

### The six fork axes, each naming the alternative NOT taken

| axis | TAKEN | NOT TAKEN (and why) |
|---|---|---|
| **functional** | Genuine `E_Legacy` (11-term, `DEFAULT_WEIGHTS`, never fitted) and genuine `E_AMBER` (ff14SB/GBn2 single point via `s20.qb2_lib.AmberSP`, bit-exactness asserted per target before any number is read) | NOT a cheap AMBER surrogate, NOT a reweighted Legacy, NOT the distogram score as a third partitioner. A surrogate would be labelled a control or not exist (BRIEF pillar 3). |
| **basis** | POINT-CLOUD Cα-RMSD of the uniform coordinate average, medoid frame, throughout | NOT built-chain, NOT projected. `readout_projected` exists and stays off; the gap is 0.156 Å of pure operator choice. |
| **readout** | Uniform coordinate average over the partition, `m` matched across partitions | NOT the medoid, NOT probability-weighted (s23 L8 closed it), NOT argmin, NOT a sharper weight. |
| **normalisation** | **AMENDED — see D1-A below.** Within-target **RANK-standardised** score of each energy across that target's own pool. Partition thresholds are **quantiles of the disagreement statistic**, fixed a priori at the quartiles below, not tuned. | NOT raw kcal/mol (the two scales differ by orders of magnitude), NOT a global z across targets (target size drives AMBER's scale), NOT a threshold fitted on RMSD, and — after D1-A — NOT the raw moment-based z-score, which is degenerate on this energy. |
| **null** | Two: (a) a **quality-matched random partition** of the same size drawn from the same pool, put through the identical readout — the operator-space control; (b) the **shared-referent floor**, i.e. the cosine two arbitrary same-size subsets of the SAME pool show against the incumbent, which is the number the disagreement cosines must beat. | NOT a uniform-on-the-torus control (a worse measure, not an uninformative one), NOT "compare to zero". Memory: `shared-referent-floor` — two quantities measured against a common reference correlate by construction; measure that floor first. |
| **THE LABEL** | The primary reported quantity is the **bias cosine against the incumbent's top-75 average**, computed by `s24/qmatch.py`'s exact recipe (`_bias` = Kabsch onto the native then subtract; `_cos`). Cα-RMSD of each partition is reported BESIDE it, descriptively. | NOT RMSD as the primary — L1 showed a screen that reads backwards can invert a whole decision, and the coordinator's instruction is explicit that the cosine is the screen that now decides. Calling this experiment "does disagreement improve RMSD" would be the wrong label and is refused. |

### D1-A — AMENDMENT TO THE NORMALISATION AXIS, 2026-09-08, BEFORE ANY OUTCOME WAS COMPUTED

**What forced it.** The pre-registered "within-target z-score" is DEGENERATE on raw AMBER.
Genuine ff14SB/GBn2 single points on **unrelaxed retrieval windows** are clash-dominated:

    target   frac E > 1e4 kcal/mol   median E      max E      frac |z_raw| < 0.1
    1CS9            0.297            -208.3      6.5e+11            0.983
    2NDN            0.610          24,772.6      7.1e+18            0.997
    5H1H              --                --            --            0.993

A single 1e18 clash sets the standard deviation, so **98.3–99.7% of every pool lands inside
|z| < 0.1** and `z(E_AMBER)` degenerates into "which candidate has the worst steric clash"
rather than "which candidate AMBER prefers". Project memory records exactly this failure:
`pauli-spectrum-delta-spike-artefact` — *a Pauli spectrum of an unconditioned energy measures
its worst clash; condition monotonically, and 99th-percentile winsorisation is not enough.*

**The amendment.** Both energies are standardised by **RANK** within the target's own pool
(`(rankdata(E) - mean) / sd`), which is a strictly MONOTONE conditioning: it changes no
ordering, no argmin and no level set, so it cannot manufacture a preference. It is also the
currency the deployed pipeline already uses (`core.pipeline._zrank`) and the same class of fix
as s20's `AMBc` control. Under rank standardisation 6.0% of candidates sit inside |z| < 0.1 —
the value expected for a well-conditioned statistic.

**This matters, and it is not cosmetic:** the `AMBER_PREFERS` partition under the two
normalisations overlaps by only **0.56–0.68**. The choice changes which candidates are studied.

**Why this is a legitimate amendment and not a fork taken after seeing the answer.** It was
diagnosed from the **energy distribution alone** — a NATIVE-FREE quantity — before a single
bias cosine, distogram alignment or RMSD was computed on any partition. No outcome variable
was inspected. The falsifier, the reading rules and my registered pessimistic prior are all
unchanged.

**Both are reported.** Rank standardisation is the declared PRIMARY; the raw moment z-score is
carried as a declared SECONDARY so the fork stays visible and the reader can see that the
degenerate arm was run rather than quietly dropped.

### D1-B — ADDED PRIMARY: ALIGNMENT WITH THE DISTOGRAM'S OWN ERROR

**L5 WAS SUBSEQUENTLY RETRACTED (2026-09-08), and this arm's reading rule is narrowed
accordingly — see `agentD_FINDINGS.md` §3.1 and §3.3.1.** Lane E's generic floor (the mean
`Dhat` over other same-length targets, not noisy) reduced L5's beta gap to −0.0107, NOT
MEASURED: a target-blind prior reproduces 102% of beta. This arm was added at the coordinator's
instruction while L5 stood, and it is **kept** because its measurement is unaffected — but what
it can support is now only *"selection aligns the output with whatever prior it is scored
against"*, not target-specific prior transfer, because this arm and its floor share the same
target-specific `Dhat` and cannot separate the two.

For every partition, alongside the bias cosine against the incumbent:

    Dhat = grid[argmin(risk, axis=1)]        the distogram's MAP distance prediction
    eP   = Dhat - Dt                         the PRIOR's own distance error
    eC   = Dc   - Dt                         the partition's emitted distance error
    cos_distogram = <eC, eP> / (|eC| |eP|)   ;  beta = <eC, eP> / |eP|^2

Taken VERBATIM from `s24/referent.py::_stats`, the coordinator's own L5 implementation, rather
than reimplemented, so the two numbers are the same object. Point-cloud basis, `min_sep=2`,
per-target then averaged across targets, no centring of either vector.

**Reading rule, fixed in advance.** The incumbent sits at ~0.640. If the Legacy-preferred and
AMBER-preferred partitions both sit at ~0.640, both physics functionals are inheriting the same
referent, and **the functional lever closes alongside the provenance lever**. A partition that
breaks materially away from 0.640 is the sprint's most valuable number and justifies the AMBER
lane on its own. Each partition's own shared-referent floor is reported beside it either way.

### D1-C — ADDED ARM: THE PHYSICS FUNCTIONAL AS A PRE-FILTER (declared before the full panel)

**The confound this removes, found on a 2-target smoke test.** The raw partitions are not
score-selected, so comparing their alignment to the incumbent's conflates *"a different
functional"* with *"no selection at all"* — and the confound is large: raw partitions came back
at 4.66–5.64 Å against the incumbent's 3.27 Å on the smoke targets. Their lower distogram
alignment could be entirely explained by not having been selected, which would be a
meaningless finding dressed as a mechanism.

**The arm.** Inside each partition, apply the **shipped distogram score** and take the same
rung (`m = min(75, |partition|)`), then the identical uniform coordinate average. Score-selected
and size-matched by construction. This is also how a physics functional would actually be
deployed — as a **pre-filter ahead of the deployed scorer** — so it is the arm that carries
deployment meaning, and it holds the scorer, the rung, the readout and the basis fixed with
only the candidate SET changing.

**This is an ADDED arm, not a changed one.** All six registered fork axes are unmoved:
functional, basis, readout, normalisation, null and THE LABEL are exactly as registered above.
It is declared here before the full 30-target panel was run.

**Reading rule.** If the pre-filtered arms land at the incumbent's RMSD and the incumbent's
distogram alignment, then filtering by genuine physics changes nothing that the distogram score
does not already determine, and the functional lever is closed. That would be a null and will
be reported as one.

### Pre-declared partition rule (not tuned)

Per target, over that target's own candidates: `zL = z(E_Legacy)`, `zA = z(E_AMBER)`,
`d = zL - zA` (positive = AMBER likes it relatively more).

    AGREE_GOOD    zL < -0.5 AND zA < -0.5          both energies rank it well
    LEGACY_PREFERS   d < q25 of d                  Legacy strongly prefers
    AMBER_PREFERS    d > q75 of d                  AMBER strongly prefers
    STRONG_DISAGREE  |d| > q75 of |d|              either direction, large magnitude
    AGREE_BAD     zL > +0.5 AND zA > +0.5          both reject it

### Reading rules, fixed in advance

* **Primary.** Mean bias cosine of each partition's emitted average vs the incumbent's, with
  the within-source control (+0.9330) and the measured shared-referent floor beside it.
  A partition **passes** only if its cosine is below ~0.65 AND its standalone mean is
  ≤ ~3.9 Å (L2's derived spec). Both, not either.
* **Secondary.** Spearman ρ(E_Legacy, E_AMBER) within target; the fraction of the classical
  top-75 each partition captures; per-partition RMSD and its own MDE.
* A cosine in the 0.85–0.95 band with no RMSD advantage is a **NULL and closes the lever**.
* 0.7–1.3× MDE is the Type-M zone and is not a result.
* Scale: AMBER is serialised under `LOCK_AMBER`, announced on open and release. The lane is
  budgeted at a target subset sized to the lock's bounded hold, and the subset is drawn by
  `s15.seed.stable_rng` and **declared before the run**, not chosen after seeing anything.

### Scope limits, declared

* **No refinement.** s23 L11: 17 of 17 restrained/unrestrained repair settings at or worse
  than no repair, bottoming out at the no-op. AMBER here is an ENERGY MEASUREMENT only. No
  minimisation, no relaxation, `steps=-1`/`k_restraint=0` single points via `AmberSP`.
* **No benchmark.** Dev-set instrument only; the sealed 60 are never touched.
* **No fitting.** No weight, threshold or partition boundary is fitted to RMSD anywhere.

---

## D2 — AN AMBER-INFORMED **SCORE** (NOT a filter, NOT a partition). PROPOSED, NOT RUN.

**Status: pre-registered and awaiting the coordinator's go. No compute spent. LOCK_AMBER not
taken.** Written now so that if it is authorised the fork list already exists.

**Why this is the last live question in Lane D, and it is a narrow one.** D1's mechanism finding
is that `rho(Legacy, distogram) = +0.3875` while `rho(AMBER, distogram) = **-0.0270**`. AMBER is
the one available functional that is *orthogonal* to the thing that already selects, and
orthogonality is the only structural reason anyone has found to expect a different functional to
help. D1 tested AMBER as a **pre-filter** and as a **partition**; both failed, and both failed
for the same reason — they let AMBER *choose candidates*, which imports AMBER's quality problem
(4.13 Å, q = 1.520) wholesale. **D2 asks whether the orthogonal direction can be inherited
without the quality being inherited**, by letting AMBER contribute to the *score* while the
distogram continues to do the choosing.

**H_D2.** A blended score `s = s_distogram + w * conditioned(E_AMBER)`, with a single global `w`
fitted under nested leave-one-fold-out CV and the same top-75 uniform readout, beats `w = 0`.

**MY REGISTERED PRIOR: I expect this to fail, and I expect it to fail for a reason that is
already on the record.** `conf.py` fitted one global exponent under nested CV over the *same*
functional and moved the endpoint by +0.0314 with 33W/93L. D1-C showed every AMBER-flavoured
candidate set costs 0.26–0.55 Å. And a blend is a monotone-preserving perturbation of an
ordering whose top-75 is already known to be insensitive to composition (L2(d), +0.0022). The
honest expectation is a null at small `|w|` and harm at large `|w|`, with nested CV selecting
`w ≈ 0`. **If nested CV picks `w = 0`, that is the result and it closes the lane.**

### Six fork axes, alternatives NOT taken

| axis | TAKEN | NOT TAKEN |
|---|---|---|
| **functional** | `s_dist + w * z_rank(E_AMBER)`, genuine ff14SB/GBn2 single points, **rank-conditioned** per D1-A | NOT raw kcal (clash-dominated, D1-A), NOT a learned AMBER surrogate, NOT AMBER replacing the distogram, NOT a per-target `w` (D1 has no power for one and §2 of the BRIEF closed routers) |
| **basis** | point-cloud Cα-RMSD, uniform coordinate average, m=75 | NOT built-chain, NOT projected |
| **readout** | the deployed uniform average, unchanged | NOT medoid, NOT probability-weighted (s23 L8), NOT a changed rung |
| **normalisation** | one **global** `w`, nested LOFO CV over the 5 pinned folds, reported with its optimism separately | NOT fitted on the folds it is scored on, NOT tuned on dev mean |
| **null** | `w = 0` (the incumbent, exactly), **plus** a sign-flipped `w` adversarial control, **plus** a rank-permuted-AMBER control that preserves the marginal but destroys the correspondence | NOT "compare to zero", NOT an unmatched arm |
| **THE LABEL** | **dev-set mean RMSD is the primary here**, unlike D1 — because a score change acts on all 126 targets and the full instrument has the power the n=30 AMBER panel did not | NOT the bias cosine as primary; it was right for D1 at n=30 and would be evasive here |

**Cost and the reason it is not free.** AMBER single points for 126 targets x 500 candidates is
~4.4x the D1 hold, i.e. roughly 40 minutes of serialised LOCK_AMBER. **I will not take the lock
without saying so first.** If authorised, the run persists `E_AMBER` per candidate to disk so
that this artefact is never re-derived at AMBER cost again — the exact failure recorded at
`agentD_FINDINGS.md` §3.3.1.

**Falsifier.** Nested-CV `w` indistinguishable from 0, or the blended arm at or worse than
`w = 0` on the dev mean. Either closes the functional lever completely.
