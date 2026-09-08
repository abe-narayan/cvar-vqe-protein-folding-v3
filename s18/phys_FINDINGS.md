# SPRINT 18 — PHYSICS findings

**Genuine Legacy (`core.energy`, eleven components, `DEFAULT_WEIGHTS`, never fitted) and genuine
AMBER ff14SB/GBn2 (`core.amber`, OpenMM). No learned surrogate stands in for either anywhere in
this workstream.**

Pre-registration: `s18/PREREG_phys.md` (§0–§6), written before any Sprint-18 physics number
existed, including my honest prior that the headline hypothesis would fail.
Coordination: `s18/COORD_phys.md`.
Modules: `s18/phys_lib.py` (the contact term in torsion space, the normalisation, the gates, the
statistics), `s18/phys_lambda.py` (Phase 6 and its declared extension), `s18/phys_decorr.py` (the
coordinator's go/no-go), `s18/phys_down.py` (Phase 8), `s18/phys_space.py` (Phase 10),
`s18/phys_report.py` (every table). Artefacts in `s18/results/`, each with an explicit `complete`
flag and row count.

Claim labels: **ESTABLISHED · SUPPORTED · PLAUSIBLE · OPEN · REFUTED · EXACT** (a theorem, not a
discovery) · **ORACLE** (needs the native).

---

---

# 0. VERDICT — leading with what most damages a physics claim, starting with my own

All five figures below are **n = 126, target as the unit, paired fold-clustered CIs, both mandatory
controls attached**, from artefacts whose `complete` flag is set.

> ## 1. MY OWN PRE-REGISTERED HEADLINE HYPOTHESIS IS REFUTED, AND ON THE OBJECTIVE THAT MATTERS IT IS REFUTED WITH A SIGN.
>
> Sprint 17 left `leg_contact` as *"real in-band information that selects 0.136 Å worse than
> random"*, and this sprint's premise was that it might succeed in the different functional role of
> a **term in the objective**. It does not.
>
> * **The pre-registered primary test** — `lam_c = +1` against `lam_c = 0`, identical start —
>   returns **+0.108 Å [−0.166, +0.377], median +0.199, 56W/70L**. Interval spans zero, median
>   worse, losing record: the pre-registered falsifier **FIRES**.
> * **On the DEPLOYED objective** (the declared extension, after MATH corrected which object the
>   brief derives and the coordinator showed the deployed form is sound) the same term at the same
>   `lam_c` and the same sign is **+0.722 Å [+0.541, +0.893], 38W / 88L — an interval excluding
>   zero.** Not "unsupported": **harmful**.
> * And it is 1.28 Å [+1.041, +1.544] worse than **doing nothing at all**, losing on 109 of 126.
>
> **The branch is closed.** No ladder extension, no re-normalisation, no sign flip.

> ## 2. AND THE REASON WAS MEASURABLE BEFORE THE LADDER WAS READ: `leg_contact` IS ORTHOGONAL TO THE DISTOGRAM'S ERROR.
>
> The coordinator's ORACLE go/no-go, at the structure the refinement starts from — the only place
> the two gradients meet:
>
>     corr(contact's descent direction, the correction the error needs)
>         = -0.006   median -0.012   vs its MJ-shuffled null +0.007 [-0.014, +0.029]   62/60
>
> **Zero on 126 targets with the interval tightened to ±0.03.** Destroying only the *direction* of
> the distogram's residual is worth ≈1.0 Å (coordinator, five controls). **That ångström is not
> reachable through this term, because this term carries no information about that direction.**
> Whatever `leg_contact`'s +0.080 [+0.032, +0.128] of in-band information is about, it is not about
> where the distogram is wrong.

> ## 3. THE RESULT THAT DAMAGES THE PROGRAMME HARDEST IS NOT ABOUT MY HYPOTHESIS. EVERY PHYSICS FILTER, PLACED IN FRONT OF THE OPERATOR THAT EMITS THE ANSWER, LOSES TO A RANDOM GATE OF THE SAME SIZE.
>
> On identical structures, each dropping the worst 50% of the same top-75 by its own score:
>
> | filter | vs a MATCHED-RANDOM gate of the same count | W/L |
> |---|---|---|
> | AMBER single point | **+0.036 [+0.006, +0.061]** | 52/74 |
> | `leg_torsion` | **+0.041 [+0.012, +0.068]** | 49/77 |
> | `leg_contact` | **+0.053 [+0.008, +0.097]** | 56/70 |
> | Legacy (11-term) | **+0.063 [+0.008, +0.125]** | 53/73 |
>
> Five scores, every interval excluding zero, every win/loss record losing. **Halving the set costs
> +0.013 Å by itself; halving it by any physics score costs four to six times that. The cost is not
> the truncation — it is the ordering.**

> ## 4. `leg_torsion`, THE ONE POSITIVE LEGACY RESULT THIS WORKSTREAM HAD, IS CLOSED BY THIS WORKSTREAM.
>
> Its near-native recall gain **reproduces** — +0.083 [+0.019, +0.135], CI excluding zero, same
> direction as Sprint 17's +0.159 [+0.073, +0.239]. `s17/phys_FINDINGS.md` §8 left it as the
> workstream's one open positive with the caveat *"the recall gain has not been shown to convert
> into any downstream accuracy."* **It is now measured on 126 targets: it does not convert. It
> converts negatively** — +0.055 [+0.026, +0.081] worse than no filter, +0.041 [+0.012, +0.068]
> worse than random.

> ## 5. A RECORDED LAW THIS PROGRAMME USES CAUSALLY IS WRONG UNDER A SELECTIVE INTERVENTION — AND I FOUND IT BY TESTING MY OWN PREDICTION.
>
> I registered, before running Phase 8, that *Legacy's gate should HELP* because
> `d_out = 1.16·d_set_mean + 0.04·d_set_best` says the operator reads what the gate improves. The
> premise held exactly — Legacy improved the set mean (3.551 → 3.527) and destroyed the set best
> (2.306 → 2.561). **The conclusion was wrong by the sign.**
>
> | gate | law PREDICTS | OBSERVED | **miss [fold CI]** |
> |---|---|---|---|
> | `leg_torsion` | **−0.024** (better) | **+0.055** (worse) | **+0.079 [+0.060, +0.098]** |
> | `legacy` | **−0.017** (better) | **+0.076** (worse) | **+0.094 [+0.063, +0.132]** |
> | **matched random** | +0.008 | +0.013 | **+0.006 [−0.003, +0.014]** |
>
> **The law holds for a random gate and fails for every score-based gate, getting the sign wrong on
> both Legacy arms.** It is a correlational fit across observed sets; it does not survive a
> *selective* intervention. Something score-ordering does to a candidate set is invisible to its
> mean and its best.

> ## 6. AND A DEFECT IN MY OWN PRE-REGISTRATION, WHICH I REPORT RATHER THAN REPAIR.
>
> `PREREG_phys.md` matched the two terms on their pool standard deviation so that `lam_c = 1` would
> read as "one sd against one sd". Measured, the **realised share** of the objective's drop supplied
> by the contact term is **0.66 at `lam_c = +0.5` and 0.85–0.94 at `lam_c = +1`**. The ladder
> therefore **never tested a small contact correction** — it tested contact-dominated objectives.
> **I am not re-normalising**, because a scale chosen after seeing an unfavourable answer is a tuned
> parameter in a methodological costume. Filed OPEN, printed beside every arm, not quietly rerun.

> ## 7. PHASE 10: THE SPACING TAX IS A DISPLACEMENT TAX AND IT IS IRREDUCIBLE — refuting this workstream's own Sprint-17 prediction.
>
> `s17/phys_FINDINGS.md` §8 item 4 predicted that repairing the coordinate average's 22.4% Cα–Cα
> contraction before AMBER would be *cheap*. Four native-free corrections, n = 126: the best,
> minimum-displacement projection onto the 3.80 Å constraint set, restores the spacing **exactly**
> at **less** Cα displacement than the incumbent (0.801 vs 0.840) and still costs **+0.023
> [+0.006, +0.037] more**. **Why**: a random displacement of the projection's own magnitude costs
> ≈+0.154 Å and the projection costs +0.163 — a null. **The bill is for moving the Cα ~0.8 Å at all,
> not for choosing a direction.** Prediction **REFUTED**, branch closed, and the AMBER pass was
> never spent.

---

## CLAIM TABLE

| claim | tier | evidence |
|---|---|---|
| `leg_contact` succeeds as a **term in the objective** where it failed as a selector | **REFUTED — my own pre-registered hypothesis** (§1.1, §3) | primary +0.108 [−0.166, +0.377] 56W/70L; on the **deployed** objective **+0.722 [+0.541, +0.893] 38W/88L** |
| the MJ table's **sequence content** matters in that role | **REFUTED, three independent measurements** (§1.2, §3) | vs residue-label permutation: −0.116 [−0.455, +0.164], −0.033 [−0.151, +0.117], +0.161 [−0.073, +0.428] — all span zero, all near-even W/L |
| a **negative** `lam_c` would rescue it | **REFUTED** (§1.2) | `lam_c = −1`: +1.314 [+0.958, +1.728], 27W/99L. Pre-declared a control, never a rescue |
| the combined objective inherits `leg_contact`'s **information** | **REFUTED — it inherits the anti-ranking** (§1.3) | ρ(objective, RMSD) monotone 0.351 → 0.282 → 0.124 → **−0.060** in `lam_c` |
| `leg_contact` is **correlated with the distogram's error** (would amplify) | **REFUTED** — it is orthogonal to it | **ORACLE**: corr = −0.006, median −0.012, null-difference +0.007 [−0.014, +0.029], n = 126 |
| the contact term raises **ensemble diversity** | **SUPPORTED but not a finding** (§1.4) | 1.053 → 2.372 — and the **zero-information control reaches 2.452**. It is the extra pair term, not the amino acids |
| refinement toward these objectives **finds** near-natives | **REFUTED** (§1.4) | sub-2 Å members per 8-start ensemble: unrefined **1.76** → full 1.21 → degree-1 0.33 → degree-1+contact **0.13** |
| **any** physics filter improves the emitted answer | **REFUTED, five scores** (§5.1) | vs no filter: +0.050 to +0.076, all CIs excluding zero, all W/L losing |
| …or beats a **matched-random gate of the same count** | **REFUTED, five scores** (§5.1) | +0.036 / +0.041 / +0.053 / +0.063 / +0.125, all CIs excluding zero |
| **Legacy's gate helps through an averaging operator** | **REFUTED — my own registered prediction** (§5.2) | premise held (set mean −0.024, set best +0.255); output **+0.076 [+0.021, +0.138]** worse |
| `d_out = 1.16·d_set_mean + 0.04·d_set_best` predicts a **gate intervention** | **REFUTED** (§5.2) | sign wrong on both Legacy arms; miss +0.079 [+0.060, +0.098] and +0.094 [+0.063, +0.132]; **holds for a random gate** (+0.006 [−0.003, +0.014]) |
| `leg_torsion`'s recall gain **converts** to downstream accuracy | **REFUTED — closing this workstream's own open item** (§5.3) | recall +0.083 [+0.019, +0.135] reproduces; output +0.055 [+0.026, +0.081] worse, +0.041 [+0.012, +0.068] worse than random |
| AMBER k = 30 is a genuine **stereochemical repair** operator | **ESTABLISHED, reproduced a fifth time** (§5.4) | vs its own gated input: clash<2 Å **5.39 → 0.00**, <2.6 Å 15.28 → 0.03, minHeavy 2.116 → **2.809**, ramaFav 0.837 → 0.875, at +0.133 Å |
| …and its **accuracy** cost | **unchanged, exactly** (§5.4) | **+0.133 [+0.112, +0.165], 24W/99L** against Sprint 17's +0.1333 [+0.1031, +0.1645], **24W/99L** |
| the AMBER convergence gate excludes the **same three targets** | **ESTABLISHED, 5th confirmation** (§5.4) | `1D6X 2NB7 7BX2`, Sprints 15/16/17×2/18 |
| AMBER repair **reorders** the filters | **REFUTED** (§5.4) | `legacy`→AMBER stays +0.078 [+0.016, +0.145] worse than `none`→AMBER; against random→AMBER the gaps shrink to null but never reverse |
| a cheaper native-free spacing restoration recovers part of the **+0.164 Å** projection tax | **REFUTED — this workstream's own Sprint-17 prediction** (§4) | best arm `minproj`: spacing 3.800 at 0.801 Å displacement, **+0.023 [+0.006, +0.037] worse than `proj`**; four arms tested |
| the spacing tax is a **direction** problem | **REFUTED — it is a displacement tax** (§4.1) | a random move of `proj`'s magnitude costs ≈+0.154; `proj` costs +0.163; the two are a null (+0.009 [−0.014, +0.033]) |
| my pool-sd normalisation means what its sentence says | **REFUTED by my own diagnostic** (§ pre-registration defect) | realised contact share **0.66–0.94**; a small contact correction is **untested**, filed OPEN |
| my torsion-space `E_contact` **is** genuine Legacy's `contact` | **EXACT** (G0) | max \|ΔE\| = **1.78 × 10⁻¹⁵** vs `s16.energy_lib`, 8 targets × 64 candidates |
| my FD torsion gradient is correct | **EXACT** (G1) | max relative error **3.5 × 10⁻¹⁰** against `s15.align_lib`'s analytic gradient |
| my instrument is the coordinator's instrument | **ESTABLISHED** (G2b) | starts bit-identical (max \|Δ\| **0.00e+00 Å**); `full` 3.606 vs `objceil` α=0 3.610, mean \|Δ\| 0.0066 |
| L30 reproduces | **ESTABLISHED** (§3) | `E_full` vs the coordinate average **+0.558 [+0.400, +0.680], 31W/95L** against L30's +0.561 [+0.407, +0.712], **31W/95L** |
| (for MATH/EXPERIMENT, not a PHYSICS claim) the **residue-additive** object matches the full objective on the real instrument | **REFUTED at n = 126** (§3) | `E_res` 4.265 vs `E_full` 3.606 — **+0.658 [+0.483, +0.828], 38W/88L**, where the 19-target enumerated instrument gave −0.004 [−0.380, +0.307] |

---

## WHAT REMAINS OPEN (short, and honestly so)

1. **A correctly-scaled *small* contact perturbation is untested.** My pre-registered normalisation
   produced a ladder that was 66–94% contact-driven at every positive rung. A scale matched on the
   terms' **realised optimisation range** rather than their pool spread is a different, legitimate
   experiment — but choosing it now, having seen this answer, would be tuning. It needs its own
   pre-registration and its own sprint, and given §2 (the term is orthogonal to the error the
   objective actually makes) I would not fund it.
2. **Why score-ordering damages an averaged set beyond its mean and its best.** §5.2 shows the
   operator law fails for score gates and holds for random gates, and the diversity column
   (`legacy` 2.213 vs `rand` 2.484 at the same m) is consistent with *correlated survivor errors*.
   Measuring that directly — the covariance structure of the survivors' errors under a score gate —
   is the one genuinely new question this lane opened, and it is cheap.
3. **Phase 10 pass 2 (AMBER on the spacing arms) was never run** and is not worth running: pass 1
   left no arm cheaper than the incumbent, which is the condition the staging was built on.

## REPRODUCTION

```
python -m s18.phys_lib                                    # gates G0, G1
python -m s18.phys_decorr                                 # the ORACLE go/no-go   ~20 min
python -m s18.phys_lambda  --out lam.json                 # PHASE 6, n = 126      ~90 min
python -m s18.phys_lambda  --bases --out lam_bases.json   # the declared extension ~25 min
python -m s18.phys_down    --out down.json                # PHASE 8, n = 126      ~2 h (AMBER)
python -m s18.phys_space   --out space.json               # PHASE 10 pass 1       ~25 min
python -m s18.phys_report                                 # every table
```
`s18/run_phys_chain.sh` runs them sequentially with BLAS pinned to one thread. Every artefact
carries `complete`, a row count and a config hash; `s18.phys_lib.read_complete` refuses to read a
partial as a result.

---

## GATES — passed before any scientific number was quoted

`python -m s18.phys_lib` → `s18/results/gates_phys.json`.

| gate | what it checks | result |
|---|---|---|
| **G0** | my torsion-space `E_contact(φ,ψ)` **is** genuine Legacy's `contact` component, against `s16.energy_lib.legacy_components_of_windows` — the call every Sprint-17 Legacy number on the record was made with | **PASS**, max \|ΔE\| = **1.78 × 10⁻¹⁵** over 8 targets × 64 candidates |
| **G1** | the batched central-difference torsion gradient is correct, measured on the pure-distance objective where `s15.align_lib` supplies an **exact analytic** gradient (so the check is of the differencing, not of a second hand-written chain rule) | **PASS**, max relative error **3.5 × 10⁻¹⁰** over 6 targets |
| **G2** | the start reproduces the record: the coordinate average and its projection | see §1 |
| **G2b** | my `full` arm and the coordinator's `alpha = 0` arm are the same objective, the same `dhat`, the same start, optimised by different code | see §1 |
| **G3** | one definition of the degree-1 object — imported from the shared adapter, never re-implemented, with the source and config hash recorded in every artefact, and a **run-time pin that refuses to continue if the definition changes mid-run** | see §1 |

**Why the gradient is finite-difference and not analytic, recorded as a choice.** The contact
energy is a function of CB, and CB is a function of (N, CA, C) of the same residue;
`core.project._torsion_grad` chains a gradient defined on CA only. Writing the N/C chain rule by
hand is exactly the kind of thing that produces a wrong number that looks right. At n ≤ 16 the
whole central-difference stencil is 4n+1 ≤ 65 structures and `core.geometry.build_backbone_batch`
builds all of them in one vectorised pass, so the analytic route buys nothing and risks a silent
defect. G1 measures the differencing rather than assuming it.

---

## THE INSTRUMENT — what "identical starting candidates" means here, stated once

Every arm in every phase starts from the same objects, so no comparison is between different sets:

* **the candidate pool** — the shipped **K = 500 BLOSUM62 retrieval pool in retrieval order**
  (`s12.instrument.load_univ`). It is **ranker-neutral**: fixed by the retrieval stage, not by any
  score being compared, so no ranker is handed a set built by its own criterion. This is the same
  set Sprint 17's 63,000 AMBER single points were computed on.
* **the filter input** — the shipped **top-75** windows (`s14.avgspace.top75_windows`), the set the
  deployed averaging operator consumes.
* **the refinement start** — the ideal-geometry projection of the coordinate average of those 75
  (`I.coordinate_average` → `I.project`). This is bit-for-bit the start `s17/refine.py` and the
  coordinator's `s18/objceil.py` use, which is what makes §1's cross-instrument check possible.
* **the distances** — `s15.distcal.gather` with the leave-fold-out `sep` bias correction, clipped
  at 2.0 Å, identical to `s17/refine.py`, `s18/objceil.py` and `s18/exp_lambda.py`.

**The two mandatory controls, and what each is for.**

| control | what it holds fixed | what it isolates |
|---|---|---|
| **zero-information** — the MJ table rebuilt on a `stable_rng` **permutation of the residue labels** | the functional form, the pair set, the switch, the magnitude distribution | the **sequence information** in `leg_contact` |
| **matched-random** — a random torsion move of the arm's own realised magnitude; or, for a gate, a random gate of the **same count** | the size of the operation | the **direction** of the move / the **ordering** the filter imposes |

A shuffled *distogram* is the right zero-information control for the **distance** term and says
nothing about the contact term; that is why the contact arms carry a shuffled **MJ table** instead.
And for a gate the matched-random control at the same count is what separates *this filter's
ordering* from the trivial effect of *averaging fewer members* — which the ledger's operator law
(`d_out = 1.16·d_set_mean + 0.04·d_set_best`) says is not neutral.

---

## AMBER DISCIPLINE, declared before use and applied everywhere below

* **The convergence gate** is *final energy finite and ≤ 1000 kcal/mol* (`core.amber.CONVERGE_MAX_KCAL`),
  declared in `PREREG_phys.md` §3 before the first minimisation ran. Every gated arm is reported
  with its **exclusion count and the excluded PDB IDs**.
* **A gated arm is compared to its own gated input.** `s18/phys_down.py` records each arm's input
  structure with the arm, and `s18/phys_lib.gated` returns the gated comparison, its ungated twin
  and the exclusion count as one inseparable object — a caller cannot get the first without the
  other two. This is the failure Sprint 17 caught in its own draft (ledger L27), where quoting a
  gated arm against the ungated instrument-wide input would have overstated a steric repair
  tenfold.
* **The protocol is the incumbent's**: harmonic restraints on N/CA/C at k = 30 kcal/mol/Å²,
  `steps = 0`, `tolerance = 1.0`, one thread, memoisation off. Sprint 17 established this set beats
  all nine alternatives at matched displacement, so it is used as found and not re-tuned.
* **Validity and RMSD are printed together, on the same structures.** `s18/phys_report.py` has no
  code path that emits a validity row without the RMSD row for the same arm, and none that emits a
  Legacy correlation without the selection delta against a matched-random operation of the same
  count. Those two substitutions are the ones the brief forbids and they are prevented structurally
  rather than by attention.

---

## PHASE 8 — the design, stated before its numbers

**The question Sprint 17 never asked.** Every physics filter in this programme has been scored as a
RANKER or a GATE, on ranking statistics. The only question a filter's user cares about is causal:
*if you put this filter in front of the operator that actually emits the answer, does the ANSWER get
better?* The ledger says the two need not agree — the terminal averaging operator consumes the set
MEAN (`d_out = 1.16·d_set_mean + 0.04·d_set_best`, R² 0.89), and Legacy's gate *raises* near-native
recall (+0.112 [+0.023, +0.198]) while *destroying* the set best (+0.219 [+0.031, +0.444]).

> **A PREDICTION, REGISTERED BEFORE THE RUN.** Put together, those two say Legacy's gate should
> **HELP** through an averaging operator: the quantity it improves is the one the operator reads,
> and the quantity it destroys enters at a weight of 0.04. If that fails, either the operator law
> is wrong or the recall gain does not reach the set mean. Either outcome is a result.

**Arms** (the brief's list): `none` · `leg_torsion` · `leg_contact` · combined `legacy` ·
`amber_sp` (a genuine ff14SB/GBn2 single point per candidate) · `rand` (matched-random gate of the
same count, 3 `stable_rng` draws) · `helix` (zero-information: Cα-RMSD to a constant ideal
α-helix). Every filter drops the worst **f = 0.50** by its own score — pre-registered as primary,
with f = 0.25 as shape — and the survivors go through the **same** averaging operator.

**Downstream operators**: the raw all-atom coordinate average (the 3.048 Å structure), its
ideal-geometry projection (the 3.213 Å incumbent), and AMBER restrained repair at k = 30. So
`legacy → amb30` and `rand → amb30` are the brief's `Legacy→AMBER` and `matched-random
gate→AMBER` arms.

---

## PHASE 10 — the design, and why it is staged

The averaging operator contracts Cα–Cα spacing to **2.949 Å against an ideal trans 3.80 Å** (22.4%,
ledger L28), and this workstream's own Sprint-17 measurement puts
Spearman(input spacing, cis fraction at `cafix`) at **−0.949** and Spearman(spacing, convergence) at
**+0.807**. The incumbent fix — the ideal-geometry projection — restores the spacing and costs
**+0.164 Å [+0.131, +0.199], 15W/111L**.

**Arms**, all native-free, all applied to the same all-atom coordinate average: `scale_uniform`
(one parameter, dilation about the centroid to mean spacing 3.80) · `bond_norm` (every Cα–Cα
direction kept, every length set to 3.80) · `segment` (bond lengths rescaled by a windowed local
target) · `minproj` (minimum displacement onto the spacing constraint set, SHAKE-style) ·
`torsion_rebuild` (ideal-geometry rebuild of the averaged torsions) · `proj` / `proj_fit` (the
incumbent) · `helix` (zero-information) · a matched-random displacement at **each arm's own realised
magnitude**.

**Staged deliberately.** Pass 1 needs no AMBER: it measures realised spacing and Cα cost, and if no
arm restores the spacing at less Cα cost than the projection the branch closes without spending a
single minimisation. Pass 2 (AMBER k = 30) runs only on what Pass 1 leaves alive. The phase is
secondary by the brief and was not allowed to take compute from Phases 6 and 8.

**A choice inside `minproj` that is not a detail.** A Cα-only correction has to be lifted to all
atoms before AMBER can see it, and re-projecting it would defeat the point — the projection is the
expensive thing being replaced. Each residue's whole atom group is therefore **translated rigidly by
its own Cα displacement**, which preserves every intra-residue bond, angle and chirality exactly and
changes only the inter-residue spacing, the quantity under test (`s18.phys_space.lift_ca`).

**Why a validity number alone would be meaningless here**, and it is why both controls are carried:
a constant α-helix has spacing 3.80, Ramachandran 1.000 and zero clashes **by construction**, and is
~4.07 Å from the native. Restoring the spacing is trivial; restoring it for less than 0.164 Å of
accuracy is the entire question.

---

## A DEFECT IN MY OWN PRE-REGISTRATION, found by measuring the thing it was supposed to control

`PREREG_phys.md` §1.1 matches the two terms on their **standard deviation over the candidate pool**,
so that `lam_c = 1` reads as *"one pool-sd of contact energy against one pool-sd of distance
objective"*. That is a defensible, native-free, pre-declared rule and I kept it. **It does not do
what the sentence implies.**

The pool sd measures how much each term varies **across retrieved candidates**. The optimiser drives
the two terms over completely different ranges: over the K = 500 pool `E_contact` has an sd of about
one unit, but a descent moves it several units, while the distance term's descent covers a small
fraction of its own pool spread. So `s18/phys_report.py` prints, beside every combined arm, the
**realised share** — the fraction of the combined objective's total drop from the start that the
contact term supplies.

> **On the smoke that first exposed this, the realised share at the PRIMARY `lam_c = +1` was 0.91–0.94.**
> A nominal "one sd against one sd" is, in the thing actually minimised, a **90%-contact objective**.
>
> **Consequence, stated plainly against my own design: the pre-registered ladder never tested a
> SMALL contact correction.** Its smallest positive rung is already contact-dominated, so what it
> measures is *"replace the objective with the contact term, keeping the distance term as a
> regulariser"*, not *"add a contact correction"*.

**And I am not fixing it by re-normalising.** A scale chosen after seeing that the first one gave an
unfavourable answer is a tuned parameter wearing a methodological costume, and this programme's
record has that failure mode written into it several times over. The ladder stands as pre-registered,
the realised share is printed beside every arm so no reader can mistake what was tested, and
"a correctly-scaled small contact perturbation is untested" is recorded as an **OPEN** question
rather than quietly run.

---

# 1. PHASE 6 — `leg_contact` AS A TERM IN THE OBJECTIVE

`s18/phys_lambda.py` → `s18/results/lam.json`, **n = 126, `complete` = true**, config hash
`ae936e1ecaad4e10`, degree-1 source `s18.exp_anova(provisional)` cfg `286331752002a94d`, pinned at
driver start and re-checked at every target (§1.0).

## 1.0 The gates, before any result

**G2 — the start reproduces the record.** Coordinate average **3.048 Å**, its ideal-geometry
projection **3.213 Å**, n = 126. Both match L30 exactly.

**G2b — cross-instrument check, and it is a strong one.** My `full` arm and the coordinator's
`alpha = 0` arm are the same objective, the same `dhat` (same leave-fold-out `sep` debias), the same
start, and different optimiser code:

| | mine | `s18/objceil.py` | agreement |
|---|---|---|---|
| **start (projection)** | 3.2126 | 3.2126 | **max \|Δ\| = 0.00 × 10⁰ Å** — bit-identical |
| **refined, α = 0 / `full`** | **3.606** | **3.610** | mean \|Δ\| **0.0066 Å**, max 0.628 on one target |

The starts are identical to the last bit and the optima agree to 0.007 Å on average, the residual
being two L-BFGS runs finding different local minima of a nonconvex objective on one target. **The
two instruments are the same instrument.**

**G3 — one definition, and a pin that enforces it.** MATH landed `s18/math_*.py` modules *while this
run was in progress*. `s18.exp_anova.build` switches to MATH's implementation the instant one exposes
the interface, which would have put the provisional `E_le1` in the first rows of the artefact and
MATH's in the rest — two definitions inside one file. `_pin_check()` records the source and config
hash at driver start and **raises rather than continue if either changes**. It did not fire; the
artefact is homogeneous. *(This hazard was not in my pre-registration. I added the guard when I saw
`math_anova.py` appear on disk mid-run, and I record that it was added mid-run rather than planned.)*

## 1.1 THE PRE-REGISTERED PRIMARY TEST — **THE FALSIFIER FIRES**

`lam_c = +1`, the single pre-declared primary arm, against `lam_c = 0` (degree-1 distance only),
identical start, n = 126, target as the unit, fold-clustered:

> ### **+0.108 Å [−0.166, +0.377], median +0.199, 56W / 70L**
>
> The interval **includes zero**, the median is **positive** (worse), and the win/loss record is
> **losing**. The pre-registered falsifier — *`lam_c = +1` does not beat `lam_c = 0` with a
> fold-clustered CI excluding zero* — **FIRES**.
>
> **H-C is REFUTED. `leg_contact` does not succeed as a term in the objective any more than it
> succeeded as a selector.** The branch closes: no ladder extension, no re-normalisation, no sign
> flip sold as a rescue.

**This is the outcome I pre-registered that I expected** (`PREREG_phys.md` §0), and the reason I
gave then is the reason the number came out this way: a weak monotone *bulk ordering* signal does not
become an argmin, and an objective's optimum is an argmin.

## 1.2 The full ladder, and the two controls

Baseline is `lam_c = 0`. `share` = the fraction of the combined objective's total drop from the start
that the contact term supplies — see the pre-registration-defect section above; it is what makes
`lam_c` interpretable.

| arm | RMSD | median | vs `lam_c = 0` [95% fold CI] | W/L | disp | **share** |
|---|---|---|---|---|---|---|
| *coordinate average* | **3.048** | 2.837 | *(the structure the pipeline builds)* | | | |
| *start (projection)* | **3.213** | 2.966 | −1.341 [−1.463, −1.154] | 100/26 | — | — |
| `lam_c = −1` **[control]** | 5.868 | 5.731 | **+1.314 [+0.958, +1.728]** | 27/99 | 5.67 | 0.94 |
| `lam_c = −0.5` **[control]** | 5.623 | 5.492 | **+1.069 [+0.713, +1.361]** | 34/92 | 5.31 | 0.88 |
| **`lam_c = 0`** | **4.554** | 4.310 | — | — | 4.62 | — |
| `lam_c = +0.5` | 4.603 | 4.319 | +0.049 [−0.105, +0.208] | 55/71 | 5.56 | 0.66 |
| **`lam_c = +1`  ← PRIMARY** | **4.662** | 4.276 | **+0.108 [−0.166, +0.377]** | **56/70** | 5.69 | 0.85 |
| `lam_c = +2` | 4.773 | 4.519 | +0.219 [−0.175, +0.597] | 57/69 | 5.90 | 0.91 |
| **`full` distance** | **3.606** | 3.400 | **−0.948 [−1.160, −0.753]** | **88/38** | 5.05 | — |
| *zero-info* shuffled-MJ at `lam_c = +1` | 4.778 | 4.498 | +0.224 [−0.085, +0.468] | 47/79 | 5.86 | — |
| *matched-random* move, same magnitude | 4.830 | 4.623 | +0.276 [+0.083, +0.494] | 46/80 | 5.56 | — |

**KEY CONTRASTS**, each as its own paired test (negative = first arm better):

| contrast | arm | base | mean [fold CI] | median | W/L |
|---|---|---|---|---|---|
| **`lam_c=+1` vs `lam_c=0`** — *the pre-registered test* | 4.662 | 4.554 | **+0.108 [−0.166, +0.377]** | +0.199 | 56/70 |
| `lam_c=+1` vs its **zero-information** control (shuffled MJ) | 4.662 | 4.778 | −0.116 [−0.455, +0.164] | −0.041 | 67/59 |
| `lam_c=+1` vs a **matched-random** move of its own magnitude | 4.662 | 4.830 | −0.167 [−0.390, +0.038] | −0.185 | 72/54 |
| **`full` vs `lam_c=0`** | 3.606 | 4.554 | **−0.948 [−1.160, −0.753]** | −0.668 | 88/38 |
| **`full` vs the START** — L30's question, reproduced | 3.606 | 3.213 | **+0.394 [+0.255, +0.523]** | +0.270 | 40/86 |
| `lam_c=0` vs the START | 4.554 | 3.213 | **+1.341 [+1.154, +1.463]** | +1.182 | 26/100 |

**Six readings, the first three of which damage my own hypothesis.**

1. **The primary test fails, as stated in §1.1.**
2. **The MJ table's SEQUENCE information is worth nothing measurable.** Against the zero-information
   control — the identical objective with the MJ table rebuilt on a permutation of the residue
   labels, so the functional form, the pair set, the switch and the magnitude distribution are all
   preserved — the primary arm buys **−0.116 Å [−0.455, +0.164], 67W/59L**. The interval includes
   zero and the win/loss is near-even. **What `leg_contact` contributes here is the *shape* of a
   pair term, not the amino-acid identities in it.**
3. **The pre-registered SIGN was right, and no negative-λ rescue was ever available.** `lam_c < 0` —
   the sign the *global* correlation ρ = −0.176 would have wanted — is catastrophic: **+1.314
   [+0.958, +1.728], 27W/99L** at `lam_c = −1`. Registering the physical sign before the run cost
   nothing and would have made a post-hoc sign flip transparently indefensible; here it is simply
   confirmed.
4. **My degree-1 base is much worse than the full objective, on my instrument, independently.**
   `full` beats `lam_c = 0` by **−0.948 [−1.160, −0.753], 88W/38L**. This is an independent
   corroboration of the coordinator's closure of the degree-1 branch, on a different module, with
   the contact machinery attached. **The base of my pre-registered primary arm is a closed object,
   and I say so rather than quoting the arm as though it were live.**
5. **L30 reproduces.** `full` is **+0.394 [+0.255, +0.523], 40W/86L** worse than the projected start
   (and 3.606 − 3.048 = +0.558 against the raw coordinate average, against L30's +0.561). *The
   pipeline still works because it does not optimise its own scoring function*, and nothing in the
   λ ladder changes that: **every** refined arm in the table is worse than the structure the
   pipeline already builds.
6. **The one thing on the correct side of the line is not a result.** The primary arm beats a
   matched-random move of its own magnitude by −0.167 [−0.390, +0.038], 72W/54L — CI includes zero.
   It is a *directional move rather than a random one*, which is the least the term could do; it is
   not a gain.

## 1.3 Objective ALIGNMENT — measured separately, and it goes the wrong way monotonically

K = 500 ranker-neutral pool, **no optimisation**. `sel` is `argmin_tied` (the tie set is averaged,
never read in cache order — the recorded trap that once invented a 1.386 Å winner).

| objective | ρ global | ρ in-band | sel global Å | sel in-band Å | **in-band sel vs MATCHED RANDOM [fold CI]** | W/L |
|---|---|---|---|---|---|---|
| `full` | **0.556** | **0.128** | 3.510 | 2.643 | +0.053 [−0.037, +0.157] | 56/67 |
| `lam_c = 0` (degree-1) | 0.351 | 0.063 | 3.974 | 2.703 | +0.113 [+0.005, +0.222] | 53/69 |
| `lam_c = +0.5` | 0.282 | 0.029 | 4.127 | 2.703 | +0.113 [+0.025, +0.209] | 52/72 |
| **`lam_c = +1`** | **0.124** | 0.039 | 4.530 | 2.741 | **+0.151 [+0.075, +0.229]** | 44/79 |
| `lam_c = +2` | **−0.060** | 0.051 | 4.906 | 2.715 | +0.125 [+0.049, +0.195] | 48/75 |
| `leg_contact` alone | **−0.176** | **0.071** | 5.634 | 2.719 | +0.130 [+0.016, +0.234] | 51/71 |
| *zero-info* shuffled MJ at `lam_c=+1` | 0.093 | −0.012 | 4.635 | 2.810 | +0.221 [+0.170, +0.280] | 38/82 |
| *matched random* | −0.000 | 0.009 | 4.235 | 2.590 | 0 by construction | — |

**Sprint 17 reproduces to three decimals on a rebuilt instrument.** `leg_contact` alone returns
ρ global **−0.176**, ρ in-band **+0.071**, global argmin **5.634 Å** — the *identical* three numbers
in `s17/phys_FINDINGS.md` §3. That is a fourth independent confirmation of that instrument.

> **Adding `leg_contact` to a distance objective destroys its global ordering monotonically in
> `lam_c`** — 0.351 → 0.282 → 0.124 → **−0.060** — while its in-band ordering stays flat and
> insignificant. **The combined objective inherits `leg_contact`'s anti-ranking, not its
> information.** The one term Sprint 17 found to carry in-band information the distance model lacks
> carries it in a form that survives neither composition nor optimisation.

## 1.4 Candidate GENERATION, near-native COVERAGE and ensemble DIVERSITY

m = 8 starts, the first top-75 windows in retrieval order — native-free and identical across arms.
`ens` = the RMSD of the ensemble's coordinate average, i.e. what a pipeline built on that objective
would emit.

| arm | set mean | set best | n < 2 Å (of 8) | **diversity** | **ens RMSD** | median | ens vs `lam_c=0` | W/L |
|---|---|---|---|---|---|---|---|---|
| **the unrefined starts** | **3.481** | **2.904** | **1.76** | 1.997 | **3.177** | 2.935 | −1.251 [−1.337, −1.146] | 102/24 |
| `lam_c = 0` | 4.558 | 4.226 | 0.33 | 1.053 | 4.427 | 4.124 | — | — |
| **`lam_c = +1`** | 4.782 | **4.013** | **0.13** | **2.372** | 4.468 | 4.235 | +0.041 [−0.339, +0.270] | 64/62 |
| `full` | 3.654 | 3.428 | 1.21 | 1.246 | 3.508 | 3.300 | −0.919 [−1.101, −0.722] | 91/35 |
| *zero-info* shuffled MJ | 4.797 | 3.973 | 0.07 | **2.452** | 4.489 | 4.182 | +0.062 [−0.163, +0.182] | 56/70 |

**Three readings, and the middle one is the only place the contact term does anything visible.**

1. **Refinement toward any of these objectives destroys near-native coverage.** The eight unrefined
   retrieval windows contain **1.76** sub-2 Å members on average; after refining toward the full
   objective **1.21**, toward degree-1 **0.33**, toward degree-1 + contact **0.13**. The pool's
   near-natives are not found by these objectives — they are *destroyed* by them.
2. **The contact term more than doubles ensemble diversity (1.053 → 2.372) and improves the set
   BEST (4.226 → 4.013) while worsening the set MEAN (4.558 → 4.782)** — the signature of a rougher
   landscape with more distinct basins. **But the zero-information control does it slightly harder**
   (diversity 2.452, set best 3.973). **The diversity is the extra pair term roughening the
   landscape, not the amino-acid identities in it.** A diversity result quoted without that control
   would have been a finding; with it, it is not one.
3. **And it does not convert.** Through the terminal averaging operator the ensemble RMSD is
   +0.041 [−0.339, +0.270], 64W/62L — a coin flip. Better set-best and higher diversity buy
   **nothing** downstream, which is the ledger's *"the terminal operator consumes the set MEAN, not
   the set BEST"* (`d_out = 1.16·d_set_mean + 0.04·d_set_best`) observed prospectively: the arm
   improved the term weighted 0.04 and damaged the term weighted 1.16.

---

# 2. THE COORDINATOR'S GO / NO-GO — **`leg_contact` IS DECORRELATED FROM THE DISTOGRAM'S ERROR**

`s18/phys_decorr.py` → `s18/results/decorr.json`, **n = 126, `complete` = true**.
**ORACLE DIAGNOSTIC** — it reads `d_true` to form the residual. It is a selector for nothing, no
parameter in this workstream is chosen with it, and it is labelled ORACLE in the module docstring
and in every line of its output.

**The question, from `s18/COORD_FINDING.md`.** With perfect distances the *identical* functional
form reaches 1.152 Å, and destroying only the **direction** of the distogram's residual buys
−1.001 [−1.226, −0.784]. So the distogram's error is structured and the structure is harmful, and
the sharp question about `leg_contact` is whether it pushes **against** that error (corrective),
**with** it (amplifying), or is **decorrelated** from it (independent).

**Definitions, with the sign nailed down** — because the obvious version of this statistic has a
sign trap in it:

    r_p    = dhat_p - dtrue_p                       the distogram's SIGNED error on pair p  [ORACLE]
    need_p = -r_p                                   the correction that pair needs
    push_p = -dE_contact/dd_p = MJ_p |switch'(d_p)| the contact term's DESCENT direction on d_p
    c_p    = MJ_p * switch(d_p)                     the literal "per-pair contribution"

    corr(push, need) > 0  CORRECTIVE     ~0  INDEPENDENT     < 0  AMPLIFYING

Every contact pair is also in the objective's pair list (**shared fraction 1.000**), so the
correlation runs over the whole contact pair set with nothing dropped. Instrument cross-check: the
distogram residual **RMS 3.205 Å, MAE 2.347 Å** — the coordinator's figures exactly.

| where | statistic | mean | median | vs its **MJ-shuffled null** [95% fold CI] | n_hi/n_lo |
|---|---|---|---|---|---|
| **at the refinement START** | **`corr(push, need)`** | **−0.006** | **−0.012** | +0.007 [−0.014, +0.029] | 62/60 |
| at the start | `corr(push, need)`, switch-active pairs only | 0.004 | 0.001 | −0.011 [−0.022, +0.005] | 53/61 |
| at the start | `corr(c, r)` — the literal statistic | 0.001 | −0.002 | −0.010 [−0.027, +0.005] | 59/63 |
| at the NATIVE | `corr(push, need)` | −0.149 | −0.155 | **+0.035 [+0.015, +0.059]** | 65/50 |
| at the native | `corr(push, need)`, switch-active only | 0.033 | 0.028 | +0.042 [−0.055, +0.114] | 56/57 |
| at the native | `corr(c, r)` — the literal statistic | **+0.175** | +0.175 | **−0.038 [−0.057, −0.023]** | 52/65 |
| over 50 pool candidates | `corr(push, need)` | −0.001 | −0.005 | +0.011 [−0.013, +0.036] | 66/60 |

*(These rows are correlations, where higher is more corrective, so the counts are printed as
`n_hi/n_lo` — the number of targets on which the real MJ table is more / less corrective than its
shuffled null. Reading them as W/L would invert the sign of every conclusion.)*

> ## The answer: the **middle** case. `leg_contact` neither corrects nor amplifies the distogram's error — at the point where it would act, it is **exactly decorrelated from it**.
>
> **At the refinement start**, which is the only place the two gradients actually meet, the
> correlation is **−0.006 with a median of −0.012**, and the difference from its MJ-shuffled null is
> **+0.007 [−0.014, +0.029]** with a 62/60 split. Over 50 pool candidates it is **−0.001**. These are
> not small effects with wide intervals; they are **zero**, on 126 targets, with the interval
> tightened to ±0.03.
>
> **Consequence, and it is a bound rather than a refutation.** The ≈1.0 Å that the coordinator's
> `shuffled` control shows is available from *not trusting the error's direction* is **not reachable
> through this term**, because this term carries no information about that direction. Whatever
> `leg_contact`'s +0.080 [+0.032, +0.128] of in-band information is about, **it is not about where
> the distogram is wrong.** That was measured *before* the λ ladder was read, and it predicts
> exactly the λ result §1 reports.

**A correction to how the statistic should be read, which I found while building it and which
changes the sign of the answer at the native.** The literal "per-pair contribution" `c_p` and the
force `push_p` have **different supports**: `c_p` is largest where the cosine switch is *saturated*,
and the force there is **exactly zero**. Correlating the energy contribution with the residual
therefore weights pairs the term **cannot move**. At the native the two disagree in sign —
`corr(c, r) = +0.175` (which reads as "correlated with the error, so it amplifies") against
`corr(push, need) = −0.149`, while restricting to the pairs the switch can actually move gives
**+0.033**. Both are reported above; **the force statistic at the start is the one that governs an
optimum**, and it is zero.

**One genuinely negative reading, stated because it is the least flattering row in the table.** At
the *native* structure `corr(push, need)` is **−0.149**, i.e. mildly amplifying. But its
MJ-shuffled null is **more** amplifying still (the difference is +0.035 [+0.015, +0.059], CI
excluding zero), and the switch-active subset is +0.033. **The amplification at the native is a
property of the pair geometry, not of the amino-acid identities** — the MJ table makes it slightly
*less* bad than a residue-label permutation does. It is not a mechanism by which `leg_contact`
could have hurt, and it is not one by which it could have helped.

---

# 3. THE DECLARED EXTENSION — the same term on the **DEPLOYED** objective, and here the answer is unambiguous

`s18/phys_lambda.py --bases` → `s18/results/lam_bases.json`, **n = 126, `complete` = true**.

**Why this exists, and that it was declared rather than fished.** Two things landed after
`PREREG_phys.md` was written, both recorded in the module's own header before it ran:

* **MATH's correction** (`s18/MATH_to_EXP.md` §2): the object the brief's ANOVA bridge derives is
  the **residue-additive** `E_res`, not the angle-additive object the provisional adapter calls
  `E_le1`; and the 19-target 2.411 Å belongs to the strict Walsh object, which has **no continuous
  image**. So my pre-registered primary arm's *base* is a coarser member of a family, and — after
  the coordinator's closure — a dead branch.
* **The coordinator's α ladder**: the deployed objective's **functional form is sound** (α = 1
  reaches 1.152 Å). That makes the sharper form of my own hypothesis — *does `leg_contact` help
  THE DEPLOYED OBJECTIVE* — the one worth answering.

**Nothing else moved.** Same contact term at `DEFAULT_WEIGHTS["contact"] = 1.0`, same pre-registered
sign, same native-free pool-sd normalisation recomputed against each base's own pool sd, same
zero-information control, same identical start. The ladder is cut to `{0, +1}` precisely so this
cannot become a fishing expedition.

| base | arm | RMSD | median | vs `lam_c = 0` [95% fold CI] | W/L | disp | share |
|---|---|---|---|---|---|---|---|
| **`E_full`** (the DEPLOYED objective) | `lam_c = 0` | **3.606** | 3.400 | — | — | 5.05 | — |
| | **`lam_c = +1`** | **4.329** | 4.136 | **+0.722 [+0.541, +0.893]** | **38/88** | 6.86 | 0.94 |
| | *zero-info* shuffled MJ | 4.361 | 4.073 | +0.755 [+0.620, +0.854] | 39/87 | 6.87 | 0.85 |
| **`E_res`** (MATH's residue-additive ANOVA) | `lam_c = 0` | 4.265 | 3.931 | — | — | 4.62 | — |
| | **`lam_c = +1`** | 4.808 | 4.515 | **+0.543 [+0.125, +0.835]** | 48/78 | 6.02 | 0.91 |
| | *zero-info* shuffled MJ | 4.647 | 4.211 | +0.382 [+0.099, +0.595] | 46/80 | 5.93 | 0.80 |

**Head-to-head contrasts:**

| contrast | | | mean [fold CI] | median | W/L |
|---|---|---|---|---|---|
| **`E_full` + contact vs `E_full`** | 4.329 | 3.606 | **+0.722 [+0.541, +0.893]** | +0.702 | 38/88 |
| `E_full` + **real** MJ vs `E_full` + **shuffled** MJ | 4.329 | 4.361 | −0.033 [−0.151, +0.117] | −0.036 | 66/60 |
| `E_res` + **real** MJ vs `E_res` + **shuffled** MJ | 4.808 | 4.647 | +0.161 [−0.073, +0.428] | +0.001 | 62/64 |
| **`E_full` + contact vs the coordinate average** | 4.329 | 3.048 | **+1.280 [+1.041, +1.544]** | +1.243 | **17/109** |
| *`E_full` alone vs the coordinate average* | 3.606 | 3.048 | **+0.558 [+0.400, +0.680]** | +0.443 | **31/95** |

> ## On the objective that is actually deployed, adding `leg_contact` is **significantly harmful**: +0.722 Å [+0.541, +0.893], 38W / 88L, interval excluding zero.
>
> The pre-registered arm's +0.108 had an interval spanning zero and could only be called "not
> supported". **On the live objective there is no ambiguity.** And the same intervention on MATH's
> residue-additive object is harmful too (+0.543 [+0.125, +0.835]), so this is not a property of one
> base.

**Three readings.**

1. **The sequence information in the MJ table is worth NOTHING, measured twice more.** Real MJ
   against a residue-label permutation is **−0.033 [−0.151, +0.117], 66W/60L** on `E_full` and
   **+0.161 [−0.073, +0.428], 62W/64L** (sign unfavourable) on `E_res`. Both span zero with
   near-even win/loss. Together with §1.2's −0.116 [−0.455, +0.164], that is **three independent
   measurements at n = 126 saying the same thing**: whatever adding a Miyazawa–Jernigan term does to
   a distance objective, it does not do it with the amino acids.
2. **L30 reproduces at the level of the win/loss record.** `E_full` alone against the coordinate
   average is **+0.558 [+0.400, +0.680], 31W/95L**, against L30's **+0.561 [+0.407, +0.712],
   31W/95L**. Same mean to 0.003 Å, same interval, and the *identical* 31/95 split. The instrument
   is the instrument.
3. **And the combined objective is 1.28 Å worse than doing nothing.** `E_full` + contact loses to
   the coordinate average on **109 of 126 targets**. Refining toward the deployed objective already
   costs 0.56 Å; adding the contact term more than doubles the damage.

**One number for MATH and EXPERIMENT, not a PHYSICS claim.** On the same start with the same
optimiser at n = 126, **`E_res` returns 4.265 Å against `E_full`'s 3.606 — +0.658 [+0.483, +0.828],
38W/88L.** On the enumerated 19-target instrument the residue-additive object was −0.004
[−0.380, +0.307] from the full objective; on the real continuous instrument at n = 126 it is
**significantly worse**. That is a fifth independent direction on the degree-1 closure, arriving
from a module built for a different purpose.

---

# 4. PHASE 10 — THE SPACING TAX IS A **DISPLACEMENT** TAX, AND IT IS IRREDUCIBLE

`s18/phys_space.py` → `s18/results/space.json`, **n = 126, `complete` = true**. Pass 1 only; Pass 2
(AMBER) was **not run**, because Pass 1 closed the branch and the staged design exists precisely so
that a closed branch costs no minimisations.

**Instrument check first.** The incumbent ideal-geometry projection costs **+0.163 [+0.132, +0.205],
24W/102L** against the coordinate average, reproducing L30's **+0.164 [+0.131, +0.199]**. The
average's Cα–Cα spacing is **2.949 Å**, reproducing L28's 2.949 against the ideal 3.80 (a 22.4%
contraction).

| arm | RMSD | median | vs the AVERAGE [fold CI] | W/L | **spacing** | **Cα disp** | vs its own MATCHED RANDOM |
|---|---|---|---|---|---|---|---|
| `none` (the coordinate average) | **3.050** | 2.822 | — | — | **2.949** | 0.000 | — |
| `scale_uniform` (1 parameter) | 4.090 | 3.719 | +1.040 [+0.871, +1.180] | 26/100 | 3.800 | 2.194 | +0.175 [+0.015, +0.347] |
| `bond_norm` (every bond → 3.80) | 3.921 | 3.765 | +0.871 [+0.715, +1.006] | 26/100 | 3.800 | 2.004 | +0.137 [−0.020, +0.308] |
| `segment` (windowed local target) | 3.919 | 3.716 | +0.869 [+0.724, +1.002] | 28/98 | 3.791 | 1.997 | +0.078 [−0.107, +0.298] |
| **`minproj`** (minimum displacement) | **3.236** | 2.993 | **+0.186 [+0.161, +0.222]** | 9/117 | **3.800** | **0.801** | **+0.038 [+0.022, +0.050]** |
| `torsion_rebuild` (averaged torsions) | 3.835 | 3.673 | +0.785 [+0.699, +0.871] | 34/92 | 3.804 | 2.174 | **−0.154 [−0.209, −0.111]** |
| **`proj`** (the incumbent) | **3.213** | 2.966 | **+0.163 [+0.132, +0.205]** | 24/102 | 3.804 | 0.840 | +0.009 [−0.014, +0.033] |
| `proj_fit` (its λ = 0 twin) | 3.205 | 2.966 | +0.155 [+0.137, +0.184] | 27/99 | 3.804 | 0.813 | — |
| *zero-info* `helix` | 4.070 | 4.243 | +1.020 [+0.773, +1.227] | 32/94 | 3.804 | 2.551 | — |

**Against the incumbent projection — the +0.164 Å reference each arm had to beat:**

| arm | vs `proj` [fold CI] | W/L |
|---|---|---|
| *do nothing* | −0.163 [−0.205, −0.132] | 102/24 |
| **`minproj`** | **+0.023 [+0.006, +0.037]** | 55/71 |
| `proj_fit` | −0.007 [−0.024, +0.007] | 67/59 |
| `torsion_rebuild` | +0.622 [+0.560, +0.716] | 43/83 |
| `segment` / `bond_norm` / `scale_uniform` / `helix` | +0.706 / +0.709 / +0.877 / +0.857 | ≈35/91 |

> ## The pre-registered success criterion — *some arm reaches the projection's spacing repair at less than +0.164 Å of Cα cost, CI excluding zero* — is **NOT MET**. Phase 10 closes.
>
> The best cheap arm, `minproj`, restores the spacing to **exactly 3.800 Å at 0.801 Å of Cα
> displacement** — *less* displacement than the projection's 0.840 — and still costs
> **+0.023 [+0.006, +0.037]** more. The interval excludes zero in the wrong direction, and the
> magnitude sits below the instrument's 0.084 Å minimum detectable effect, so the honest reading is
> **"matches the projection, does not beat it"**, not "is worse than it".

## 4.1 WHY — and this is the part worth keeping

Read the **matched-random** column, which is what the pre-registration carried it for. A random
displacement of the *same magnitude* as the projection costs about **+0.154 Å** against the average
(the projection is +0.009 [−0.014, +0.033] away from it — a null). A random displacement of
`minproj`'s magnitude costs about **+0.148 Å** (`minproj` is +0.038 worse than it).

> **So the cost of restoring the spacing is not the cost of choosing a particular direction. It is
> the cost of MOVING THE Cα BY ~0.8 Å AT ALL.** Undoing a 22.4% contraction on a 3.80 Å backbone
> *requires* ~0.8 Å of Cα displacement, and ~0.8 Å of displacement away from the coordinate average
> costs ~0.16 Å of Cα-RMSD **whichever way you go**. The projection is not paying for being clever;
> it is paying the geometric bill, and it happens to be at the cheap end of the distribution.
>
> **This closes the fourth open item of `s17/phys_FINDINGS.md` §8** — the prediction that repairing
> the coordinate average's spacing *before* AMBER would be cheap. It is not. The prediction was this
> workstream's own, it is now measured at n = 126, and it is **REFUTED**.

**Two smaller readings, one of which is against the arms I built.**

1. **Only one arm beats its own matched random, and it is the most expensive one.**
   `torsion_rebuild` — the ideal-geometry rebuild of the *averaged torsions* — is
   **−0.154 [−0.209, −0.111]** better than a random displacement of its own 2.17 Å magnitude. It is
   genuinely a *purposeful* move, and it still lands at +0.785 Å against the average, because the
   move is three times too big. Every cheap arm I built (`scale_uniform`, `bond_norm`, `segment`,
   `minproj`) is at or **worse than** its matched random: their directions carry no information the
   magnitude does not already explain.
2. **The projection's Ramachandran penalty is free.** `proj` (λ = 0.3 Ramachandran-penalised) and
   `proj_fit` (λ = 0) are **−0.007 [−0.024, +0.007], 67W/59L** apart — indistinguishable. The
   penalty buys the validity it is there for at no measurable Cα cost.

---

# 5. PHASE 8 — THE CAUSAL DOWNSTREAM COMPARISON: **EVERY PHYSICS FILTER MAKES THE ANSWER WORSE, AND A RANDOM GATE OF THE SAME SIZE BEATS ALL OF THEM**

`s18/phys_down.py` → `s18/results/down.json`, **n = 126, `complete` = true**. Every filter sees the
**same** shipped top-75 windows and drops the worst 50% by its own score; the survivors go through
the **same** averaging operator. f = 0.50 is the pre-registered primary; f = 0.25 is shape.

## 5.1 The emitted answer

| filter | avg RMSD | median | **vs NO FILTER [fold CI]** | W/L | proj | m | set mean | set best | div |
|---|---|---|---|---|---|---|---|---|---|
| **`none`** | **3.050** | 2.822 | — | — | **3.213** | 75 | 3.551 | **2.306** | 2.486 |
| `leg_torsion` | 3.105 | 2.902 | **+0.055 [+0.026, +0.081]** | 46/80 | 3.255 | 38 | **3.523** | 2.501 | 2.214 |
| `leg_contact` | 3.116 | 2.862 | **+0.067 [+0.014, +0.119]** | 54/72 | 3.275 | 38 | 3.589 | 2.468 | 2.489 |
| `legacy` (11-term) | 3.126 | 2.954 | **+0.076 [+0.021, +0.138]** | 49/77 | 3.276 | 38 | **3.527** | 2.561 | 2.213 |
| `amber_sp` | 3.099 | 2.922 | **+0.050 [+0.011, +0.081]** | 53/73 | 3.252 | 38 | 3.553 | 2.444 | 2.394 |
| **`rand` [CONTROL]** | **3.063** | 2.855 | **+0.013 [+0.003, +0.024]** | 56/70 | 3.220 | 38 | 3.554 | 2.422 | 2.484 |
| `helix` [ZERO-INFO] | 3.188 | 3.131 | +0.139 [+0.069, +0.193] | 48/78 | 3.329 | 38 | 3.538 | 2.749 | 2.036 |

**Against the matched-random gate of the same count — the comparison that separates *this filter's
ordering* from *averaging fewer members*:**

| filter → | vs `rand` gate | median | W/L |
|---|---|---|---|
| `amber_sp` | **+0.036 [+0.006, +0.061]** | +0.012 | 52/74 |
| `leg_torsion` | **+0.041 [+0.012, +0.068]** | +0.018 | 49/77 |
| `leg_contact` | **+0.053 [+0.008, +0.097]** | +0.018 | 56/70 |
| `legacy` | **+0.063 [+0.008, +0.125]** | +0.025 | 53/73 |
| `helix` [zero-info] | +0.125 [+0.062, +0.181] | +0.057 | 44/82 |

> ## **Every physics filter in the programme, placed in front of the operator that emits the answer, is significantly WORSE than a random gate that removes the same number of candidates.** Five scores, five intervals excluding zero, five losing win/loss records.
>
> Halving the set costs **+0.013 Å** by itself. Halving it *by any physics score* costs **+0.050 to
> +0.076 Å** — four to six times as much. **The cost is not the truncation. It is the ordering.**

**And it reproduces at f = 0.25** (shape, not promoted): `rand` +0.002 [−0.001, +0.005] — free —
against `leg_torsion` +0.042, `leg_contact` +0.029, `legacy` +0.028, `amber_sp` +0.029, `helix`
+0.068. Same sign, same ranking, smaller magnitude.

## 5.2 **MY REGISTERED PREDICTION IS REFUTED — and the recorded operator law fails with it**

`s18/phys_down.py`'s header registered, before the run: *Legacy's gate should HELP through an
averaging operator, because the ledger's law `d_out = 1.16·d_set_mean + 0.04·d_set_best` says the
operator reads the quantity Legacy improves (the class, hence the set mean) and only barely reads the
quantity it destroys (the set best).*

**The premise held. The conclusion did not.** Legacy's gate *did* improve the set mean
(3.551 → 3.527) and *did* destroy the set best (2.306 → 2.561), exactly as Sprint 17 measured. The
output still got **+0.076 Å worse**.

**So I tested the law directly, as a causal predictor, arm by arm:**

| filter | Δ set mean | Δ set best | **law PREDICTS** | **OBSERVED** | **miss [fold CI]** |
|---|---|---|---|---|---|
| `leg_torsion` | −0.028 | +0.195 | **−0.024** (better) | **+0.055** (worse) | **+0.079 [+0.060, +0.098]** |
| `legacy` | −0.024 | +0.255 | **−0.017** (better) | **+0.076** (worse) | **+0.094 [+0.063, +0.132]** |
| `amber_sp` | +0.003 | +0.137 | +0.008 | +0.050 | +0.041 [+0.024, +0.062] |
| `leg_contact` | +0.038 | +0.162 | +0.051 | +0.067 | +0.016 [−0.001, +0.034] |
| **`rand` [CONTROL]** | +0.003 | +0.116 | **+0.008** | **+0.013** | **+0.006 [−0.003, +0.014]** |
| `helix` [zero-info] | −0.013 | +0.442 | +0.003 | +0.139 | +0.135 [+0.103, +0.164] |

> ### The law **holds for a random gate** (miss +0.006, CI spanning zero) and **fails for every score-based gate**, getting the *sign* wrong on the two Legacy arms.
>
> `d_out = 1.16·d_set_mean + 0.04·d_set_best` was fit across *observed* candidate sets. It does not
> survive a **selective** intervention. Something a score-ordered subset does to the coordinate
> average is invisible to its mean and its best — the natural reading is that scoring **correlates
> the survivors' errors**, and the averaging operator's accuracy comes from error *cancellation*
> (ledger L30), so a set whose errors agree averages worse than a set of the same mean quality whose
> errors are independent. The diversity column is consistent with it (`legacy` 2.213 and
> `leg_torsion` 2.214 against `rand`'s 2.484 at the same m) but does not prove it.
>
> **Recorded as a correction to a law this programme has been using as though it were causal.**

## 5.3 `leg_torsion`'s recall gain reproduces — **and does not convert**

Near-native recall, on the **55 targets that have a sub-2 Å member**, one denominator, against a
matched-random gate of the same count:

| filter | recall | random | **recall − random [fold CI]** |
|---|---|---|---|
| **`leg_torsion`** | **0.568** | 0.485 | **+0.083 [+0.019, +0.135]** |
| `amber_sp` | 0.547 | 0.485 | **+0.062 [+0.032, +0.098]** |
| `legacy` | 0.509 | 0.485 | +0.024 [−0.000, +0.054] |
| `leg_contact` | 0.476 | 0.485 | −0.009 [−0.064, +0.054] |
| `helix` [zero-info] | 0.493 | 0.485 | +0.008 [−0.072, +0.106] |

> **`leg_torsion` is still the one Legacy gate that raises near-native recall** — +0.083
> [+0.019, +0.135], CI excluding zero, the same direction as Sprint 17's +0.159 [+0.073, +0.239]
> (smaller here because this gate operates on the already-filtered top-75 rather than K = 500).
>
> **And it still makes the emitted answer +0.055 [+0.026, +0.081] worse, and +0.041 [+0.012, +0.068]
> worse than a random gate.** `s17/phys_FINDINGS.md` §8 left this as the workstream's one open
> positive — *"not a result yet: the recall gain has not been shown to convert into any downstream
> accuracy."* **It is now measured. It does not convert. It converts negatively.** The last
> surviving positive Legacy result in this programme is closed by its own workstream.

## 5.4 AMBER at k = 30 — the repair reproduces exactly, and it does not reorder the filters

**THE CONVERGENCE GATE**, declared before use (final energy finite and ≤ 1000 kcal/mol), reported
with its exclusions **by name**. Every arm is compared to **its own gated input**.

| arm | gate excl | excluded | RMSD | **its own gated input** | **vs that input [fold CI]** | W/L | Cα disp |
|---|---|---|---|---|---|---|---|
| `none` → AMBER | **3** | `1D6X 2NB7 7BX2` | 3.183 | 3.049 | **+0.133 [+0.112, +0.165]** | **24/99** | 0.713 |
| `leg_torsion` → AMBER | 3 | `1D6X 2NB7 6QAX` | 3.218 | 3.097 | +0.121 [+0.099, +0.150] | 19/104 | 0.638 |
| `leg_contact` → AMBER | 1 | `9KAR` | 3.210 | 3.086 | +0.124 [+0.108, +0.151] | 22/103 | 0.700 |
| **`legacy` → AMBER** | 1 | `8TXS` | 3.260 | 3.127 | +0.133 [+0.100, +0.173] | 19/106 | 0.653 |
| **`rand` → AMBER** [CONTROL] | 2 | `1G89 1MF6` | 3.211 | 3.070 | +0.141 [+0.107, +0.185] | 23/101 | 0.718 |

> **The instrument reproduces to the win/loss record.** The unfiltered arm returns
> **+0.133 [+0.112, +0.165], 24W/99L** against Sprint 17's **+0.1333 [+0.1031, +0.1645], 24W/99L** —
> same mean to 0.0003 Å, same interval, the *identical* 24/99 split — and the gate excludes the
> **same three targets, `1D6X 2NB7 7BX2`**, for the fifth independent time (Sprints 15, 16, 17 ×2,
> 18). This is as strong an instrument check as this programme has.

**VALIDITY, on the same structures whose RMSD is quoted, each arm above its own gated input:**

| structure | ramaFav ↑ | ramaOut ↓ | clash<2 ↓ | clash<2.6 ↓ | minHeavy ↑ | cis ↓ | ωdev ↓ |
|---|---|---|---|---|---|---|---|
| `none` → AMBER | **0.875** | 0.041 | **0.00** | **0.03** | **2.809** | 0.074 | 21.5 |
| ↳ *its own gated input* | 0.837 | 0.080 | 5.39 | 15.28 | 2.116 | 0.000 | 0.0 |
| `leg_torsion` → AMBER | 0.894 | 0.037 | 0.00 | 0.00 | 2.833 | 0.056 | 18.7 |
| ↳ *its own gated input* | 0.853 | 0.065 | 5.51 | 14.89 | 2.135 | 0.000 | 0.0 |
| `legacy` → AMBER | 0.890 | 0.045 | 0.00 | 0.01 | 2.805 | 0.051 | 17.7 |
| ↳ *its own gated input* | 0.873 | 0.059 | **6.79** | **17.81** | 2.062 | 0.000 | 0.0 |
| `rand` → AMBER [CONTROL] | 0.887 | 0.037 | 0.00 | 0.02 | 2.795 | 0.077 | 23.2 |
| ↳ *its own gated input* | 0.846 | 0.080 | 5.71 | 16.47 | 2.050 | 0.000 | 0.0 |

**AMBER does the one job it has earned, on every arm**: clashes below 2.0 Å to **0.00**, below 2.6 Å
to ≤0.03, minimum heavy separation +0.69 Å, Ramachandran-favoured *up* 0.837 → 0.875, at ~0.7 Å of
Cα displacement and +0.13 Å of Cα accuracy. **And that is a validity statement, not an accuracy one —
every arm's RMSD is printed beside it and every arm's RMSD gets worse.**

**Does repair reorder the filters? No.**

| | | | mean [fold CI] | median | W/L |
|---|---|---|---|---|---|
| `legacy`→AMBER vs `none`→AMBER (n = 122) | 3.265 | 3.188 | **+0.078 [+0.016, +0.145]** | +0.027 | 43/79 |
| `leg_torsion`→AMBER vs `none`→AMBER (n = 122) | 3.220 | 3.176 | **+0.043 [+0.012, +0.073]** | +0.005 | 58/64 |
| `legacy`→AMBER vs `rand`→AMBER (n = 123) | 3.266 | 3.214 | +0.052 [−0.018, +0.128] | +0.010 | 56/67 |
| `leg_torsion`→AMBER vs `rand`→AMBER (n = 121) | 3.225 | 3.200 | +0.025 [−0.015, +0.063] | −0.002 | 63/58 |

**Legacy→AMBER stays significantly worse than no-filter→AMBER**; against the matched-random
gate→AMBER the differences shrink to null, i.e. the repair *washes out* part of the filter's damage
but never reverses its sign. **The brief's `Legacy→AMBER` and `matched-random gate→AMBER` arms
therefore give the same verdict as everything else in this phase.**

## 5.5 Two mechanistic notes from the validity table

1. **Legacy's gate selects for a WORSE-PACKED average.** The un-repaired coordinate average of
   Legacy's survivors carries **6.76** sub-2.0 Å clashes and a minimum heavy separation of
   **2.057 Å**, against no-filter's **5.61** and **2.096** and a random gate's 5.63 and 2.051 —
   while `leg_contact`'s survivors give the *best*-packed average (4.77 clashes, 2.151 Å). An
   eleven-term potential whose largest weight is `steric = 4.0` produces, after averaging, a
   structure with **more** steric violations than no selection at all. That is the set-mean trap
   showing up in a physical observable rather than in an RMSD.
2. **Halving the set reduces the averaging contraction, slightly.** Cα–Cα spacing rises from
   **2.949 Å** (75 members) to **2.964–3.058 Å** (38 members) on every arm including random —
   fewer members, less error cancellation, less contraction. Consistent with L28's mechanism, and
   it is *not* a benefit: the arms with the most spacing recovery (`leg_torsion` 3.058, `helix`
   3.060) are also among the worst on RMSD.
