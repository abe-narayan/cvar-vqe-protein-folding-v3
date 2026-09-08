# SPRINT 16 — RETRACT findings

Agent: `retract`. Date: **2026-09-06**. Two correction tasks handed over by the LIT workstream.

Modules: `s16/retract_law.py`, `s16/retract_exact.py`, `s16/retract_screen.py`.
Artefacts: `s16/results/retract_law.json`, `retract_exact.json`, `retract_screen.json`,
`retract_exact.log`. Also written to `s12/results/s16_retract_exact.json` via `I.write`.

Tiering: **DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED**.
Statistical unit is the **TARGET**, n = 126, median printed beside mean, and every interval is
given twice — `I.paired`'s i.i.d.-target bootstrap **and** a two-stage fold-clustered bootstrap
(resample the 5 folds, then targets within), because `I.paired` does not cluster and BRIEF §3.4
asks for one that does. **With only 5 folds the clustered interval is coarse; it is reported, not
relied on.** Seeded with `s15/seed.py`'s `stable_rng`. No `hash()`. The **sealed 60-target
benchmark was not read, probed, or referenced.** Natives are used here only for **auditing a
published number** — permitted by BRIEF §3.1 — and nothing in this document enters a predictor.

---

## THE HEADLINE, worst-first

> **1. Sprint 15's "parameter-free fusion law, validated to 0.162 Å on 126 targets" predicts
> nothing. It is an exact algebraic identity, and the 0.162 Å it was "validated to" is 54% an
> arithmetic-mean error and 46% a frame-convention mismatch — with 0% attributable to the identity
> failing.** In a single common frame with the correct mean, the residual on all 126 real targets is
> **2.65 × 10⁻¹⁵ Å**.
>
> **2. LIT's two "surviving novelties" both fail.** (a) "It holds through Kabsch superposition on
> 126 real targets" is a **theorem**, not an empirical finding — and Sprint 15 never measured it,
> because it measured a different frame. (b) "`s` is native-free, so it works as an *a priori*
> fusion screen" is **REFUTED**: the law's *prediction* needs `r`, a distance to the native, and `s`
> alone ranks "does fusion win on this target" at **AUC 0.401 against a 0.500 null** — it points the
> **wrong way**.
>
> **3. The Sprint 15 numerical audit passed this row.** `s15/verify_FINDINGS.md` verified that every
> published fusion number traces to `coherence.json` **exactly** — which it does. A transcription
> audit cannot see a formula error. That is the transferable lesson, and it is recorded in the audit
> file itself.
>
> **4. Shao & Zhu (2018) is reinstated, and two of the natural defences against it fail.** It tests
> **ff14SB and FF14SBonlysc with GB-Neck2** — our exact pair. The "they only tested longer peptides"
> defence fails (**1E0Q is a 17-residue β-hairpin**, one residue above our band, and is one of the
> systems reported discrepant). The "free energy is a weaker observable than potential-energy
> ranking" defence is **backwards** — free energy is Boltzmann-weighted and therefore stronger.
> **"AMBER anti-ranks peptide natives" is documented expected behaviour and may not be presented as
> a discovery.**

---

# TASK A — THE FUSION LAW

## A.1 The inventory, complete, produced before anything was changed

**Computation (1 module, 2 artefacts).**

| file | role | status |
|---|---|---|
| `s15/coherence.py` line 123 (`r_bar = 0.5 * (rd + rp)`) | **the defect** — arithmetic mean where the identity needs the quadratic | corrected in place; both means now emitted, original retained as `law_prediction` |
| `s15/results/coherence.json` | the n = 126 artefact the whole record quotes | per-target `rows` **correct and reused**; aggregate `law_prediction` / `prediction_error` **superseded** |
| `s12/results/s15_coherence.json` | `I.write` mirror, `complete: true`, `n_rows: 126` | same |

**A provenance defect found in passing.** `s15/results/coherence.json` on disk carries the key
`coordavg_vs_best_single`, which the mislabelling correction inside `coherence.py` renamed to
`coordavg_vs_better_channel_NATIVE_FREE`. **The persisted artefact was therefore written by an
earlier version of the module than the one in the tree, and has not been regenerated.** Its rows are
sound (verified below, bit-identically), but the file is not a product of the current source. Noted
in the module's correction block.

**Documents quoting the law (9 files, 13 sites) — all now carry a dated Sprint 16 correction with
the original preserved, struck through, never deleted:**

| file | site |
|---|---|
| `s15/coord_FINDINGS.md` | §"The fusion law" (n = 126 table); §"The law" (n = 24 derivation and its unnecessary caveats); §"Why this is worth having" (the screen claim); the numbered claim list |
| `s15/FINAL_DOSSIER.md` | §1.5, and the contributions list |
| `s15/RESULTS.md` | §6.2 |
| `s15/INFORMATION_CHANNEL.md` | §4 table and the `s`-is-native-free paragraph |
| `s15/PAPER_DRAFT.md` | abstract "**Fourth**"; §6 in full; the intro summary; the figure list |
| `s15/REVIEWER_RESPONSES.md` | **R3.4** (which conceded, but defended the wrong thing) and "what the reviews changed" item 3 |
| `s15/LEGACY_VS_AMBER.md` | the fusion paragraph |
| `s15/PROTOCOL_FROZEN.md` | the "fusion law's calibration (24 targets)" exclusion |
| `s15/figures/README.md` + `s15/figures.py` | `fig12_fusion_law` — **plots `law_prediction`, i.e. the wrong mean. Marked "do not publish as drawn."** |

**Not edited, flagged for the coordinator:** `s16/BRIEF.md` §2 still lists *"fusion law
`d_avg ≈ √(r² − (s/2)²)`, parameter-free, predicts to 0.162 Å on 126 targets, and s is
native-free"* among the established Sprint 15 results. That line is now wrong on three counts. The
BRIEF is the shared contract and coordinator-owned, so it was left alone rather than edited
mid-sprint.

## A.2 The corrected numbers — **DEMONSTRATED**

`s16/retract_law.py`, n = 126, **no refit**: the law consumes only `rmsd_disto`, `rmsd_pool`,
`disagreement` and `rmsd_coordavg`, all persisted per target. The arithmetic reconstruction
reproduces the recorded `law_prediction` to `atol 1e-9` on every target (asserted in the module), so
the correction is applied to the identical inputs.

| | ~~arithmetic mean (published)~~ | **quadratic mean (correct)** |
|---|---|---|
| mean law prediction | ~~3.2049~~ | **3.2919** |
| measured coordinate average | 3.3668 (median 3.1668) | 3.3668 |
| **mean signed residual** | ~~**+0.1619**~~ | **+0.0748** |
| **median signed residual** | ~~+0.0840~~ | **+0.0186** |
| targets where the law **under**-predicts | 126 / 126 | 126 / 126 |

**Paired difference in |residual|, quadratic − arithmetic, TARGET as the unit:**

| | value |
|---|---|
| mean | **−0.0871** |
| i.i.d.-target bootstrap 95% CI | **[−0.1185, −0.0605]** |
| **fold-clustered** bootstrap 95% CI | **[−0.1285, −0.0556]** |
| median | −0.0150 |
| W/L | **126 / 0** |
| per-fold | −0.091, −0.122, −0.100, −0.074, −0.056 — **5/5 same sign** |

The direction is forced: the quadratic mean is ≥ the arithmetic mean, so the prediction can only
rise, and the law under-predicts on every target. **54% of Sprint 15's reported "prediction error"
was an arithmetic artefact.**

**LIT's diagnosis of *when* it matters was wrong, and in the direction that understates the
problem.** LIT wrote that "they agree here only because r₁ ≈ r₂" and that "this particular row is
not affected much". That is true of the **aggregate means** (3.622 vs 3.659) and false of the
**targets**, which are the unit:

| |r₁ − r₂| | value |
|---|---|
| mean / median | **0.980** / 0.567 Å |
| p75 / p90 / p95 / max | 1.694 / 2.343 / 2.673 / **4.798** Å |
| fraction > 1 Å / > 2 Å | **38.1%** / **13.5%** |

correlation of the gap with the size of the substitution: **+0.844**. Worked examples:

| target | r₁ | r₂ | gap | measured | ~~err arith~~ | **err quad** |
|---|---|---|---|---|---|---|
| 2BP4 | 5.758 | 0.960 | 4.798 | 3.200 | ~~+1.403~~ | **+0.202** |
| 6RRO | 3.574 | 0.405 | 3.169 | 1.860 | ~~+0.901~~ | **+0.008** |
| 2RUO | 4.727 | 1.203 | 3.523 | 2.787 | ~~+0.763~~ | **+0.104** |

## A.3 Is the identity exact through Kabsch superposition? — **DEMONSTRATED: yes, and it is a theorem**

This was the one part of the Sprint 15 result LIT thought might be non-trivial, because Kabsch
superposition is a minimisation and therefore nonlinear, and because `s16/STEER_FINDINGS.md` has
just shown that this programme's errors live outside the linear regime. **It is not non-trivial.**

`s16/retract_exact.py` re-runs `s15/coherence.py`'s two fits with the same starts, seeds and
optimiser. **Reproduction is bit-identical: max |Δ| = 0.000e+00 against the persisted artefact on
all four quantities** (`rmsd_disto`, `rmsd_pool`, `disagreement`, `rmsd_coordavg`). Four frame
conventions, n = 126, quadratic mean:

| arm | what it does | mean \|residual\| | max | sign |
|---|---|---|---|---|
| **A1 `S15_pipeline`** | Sprint 15's convention — five independent superpositions | 7.48e−02 | 6.26e−01 | 126/0 |
| **A2 `common_frame`** | both structures superposed onto the native once; `s` and `d_avg` are plain RMSDs in that frame | **2.65e−15** | 2.52e−14 | **55/61** |
| A3 `s_common_only` | pipeline `d`, common-frame `s` | 1.74e−01 | 1.98e+00 | 126/0 |
| A4 `d_common_only` | pipeline `s`, common-frame `d` | 9.90e−02 | 1.61e+00 | 0/126 |
| A5 `commonavg_then_kabsch` | common-frame average, then Kabsch onto the native | **4.24e−15** | 6.21e−14 | 52/64 |

**A2 is machine epsilon with a balanced sign (55 under / 61 over) — floating-point noise, not a
residual.** LIT's exactness claim is confirmed on the real data.

**And the same run isolates the whole of Sprint 15's residual as a frame convention:**

| source | size |
|---|---|
| averaging in the **medoid** frame instead of the target's (`I.coordinate_average` superposes onto the medoid) | **+0.1738** Å |
| `s` being **Kabsch-minimised** (3.138 Å) rather than measured in that same frame (3.506 Å) | **−0.0990** Å |
| = the A1 residual | **+0.0748** Å, exactly |
| the final Kabsch of the average onto the native | **0.000000** (mean = median = max, 126/126) |
| the arithmetic-mean artefact, on top | **+0.0871** Å |
| = the published figure | **+0.1619** Å |

### The lemma, which is why "it holds through Kabsch" was never in doubt

The final superposition contributes **exactly zero** because of a two-line fact. Kabsch's optimality
condition for `R = I` is that the cross-covariance `A = Σₙ Wₙ Tₙᵀ` is symmetric PSD. `A` is **linear
in the moving structure**, so for structures `W⁽¹⁾…W⁽ᴮ⁾` each already optimally superposed onto `T`,
the mean structure has `A = (1/B) Σ_b A_b` — a mean of symmetric PSD matrices, hence symmetric PSD.
**The average of optimally-superposed structures is itself optimally superposed.** Verified: 0
violations in 3,000 random trials (`stable_rng`, random `B ∈ [2,6)`, `n ∈ [6,20)`) and exactly
0.000000 on all 126 real targets.

> **Therefore LIT's surviving novelty (a) — "it holds through Kabsch superposition on 126 real
> targets" — is not a finding. It is forced, and Sprint 15 did not measure it.** What Sprint 15
> measured, and reported as a prediction error, was the difference between two frame conventions
> plus the wrong mean.

**ORACLE DIAGNOSTIC, recorded because it is interesting and unusable.** Averaging in the *native's*
frame instead of the medoid's gives **3.1930 Å against 3.3668 Å — 0.174 Å better**. It requires the
native. It is a ceiling, not a method, and is labelled as such wherever it appears.

## A.4 Is `s` an a-priori native-free fusion screen? — **REFUTED**

LIT's surviving novelty (b), tested rather than repeated (`s16/retract_screen.py`, n = 126, no
refit).

**First, the claim as written is not available at all.** The law predicts `d_avg` from `s` **and
`r`**, and `r` is a mean RMSD **to the native**. Only the *input* `s` is native-free; the
*prediction* is not. So "the fusion decision becomes a calculation on unlabelled data" is false on
its face.

**Second, the mechanism is wrong.** The ambiguity decomposition guarantees the ensemble beats the
**average** member by the diversity, for *any* `s > 0`. It is silent on beating the **best** member —
which is the only question a screen is for. **The identity answers a question nobody asks.**

**Third, measured:**

| | value |
|---|---|
| corr(`s`, gain over the **quadratic-mean member**) | +0.777 (Spearman +0.770) — **algebra, not prediction**: the identity makes this gain a deterministic function of `s` and `r` |
| corr(`s`, gain over the **better single channel**) — the practitioner's quantity | **+0.112** (Spearman **+0.040**) |
| **AUC of `s` as a ranker of "fusion wins on this target"** | **0.401** against a 0.500 null — **anti-informative** |
| mean `s` where fusion wins / loses | 2.947 / **3.493** — higher disagreement predicts fusion *losing* |
| leave-fold-out threshold "fuse iff `s > t`", `t` chosen **supervised** on 4 folds | every fold chose the grid minimum, i.e. **"always fuse"**, 5/5 |
| screen − always-fuse | +0.0002 [+0.0000, +0.0006] i.i.d.; [+0.0000, +0.0009] fold-clustered |
| always-fuse − never-fuse | −0.2552 [−0.3690, −0.1442] i.i.d.; [−0.4385, −0.0986] fold-clustered, median −0.1384, W/L 81/44 |
| ORACLE per-target pick (a ceiling) | 3.2420 vs always-fuse 3.3668 — **0.125 Å of headroom exists; `s` does not reach it** |

The threshold is **fold-honest supervised, not native-free** — it minimises held-out-fold RMSD,
exactly the status `s16/LEDGER.md` L8 assigns to the steering step size — and it is labelled as such.
Since "always fuse" is on the grid, a training fold can always decline to screen, and it always did.

**Fusion itself is worth −0.255 Å over the better single channel and that is unaffected.** What is
refuted is the *screen*.

## A.5 What survives of the Sprint 15 fusion result

| claim | verdict |
|---|---|
| "a parameter-free law" | **WITHDRAWN** — Krogh & Vedelsby (1995) ambiguity decomposition, two-member case. Cite it |
| "predicting an unseen quantity to 0.162 Å across 126 targets" | **WITHDRAWN** — it predicts nothing unseen; the residual is 0.075 Å with the right mean and 2.65e−15 in the right frame |
| "exact when the native is equidistant and the three points are coplanar" | **WITHDRAWN** — both caveats unnecessary; exact in any inner-product space |
| "it holds through Kabsch superposition on 126 real targets" (LIT's (a)) | **TRUE BUT NOT A FINDING** — a theorem, and Sprint 15 measured a different frame |
| "`s` is native-free, hence an a-priori fusion screen" (LIT's (b)) | **REFUTED** — AUC 0.401 vs a 0.500 null; the prediction needs `r`, which is not native-free |
| per-pair error correlation **+0.602** coexisting with **3.138 Å** structural disagreement | **SURVIVES — DEMONSTRATED.** A measurement, not an identity. This is the durable content |
| coordinate averaging beats the better single channel by −0.255 [−0.370, −0.145] | **SURVIVES** (unexamined here; unaffected by the correction) |
| the `s/r` scaling statement ("an ångström needs `s → 2r`") | **SURVIVES** — correct algebra, and now correctly attributed |

**Net novelty of the fusion result: none.** What remains is one measurement (+0.602 vs 3.138 Å) and
one negative result that is genuinely worth reporting — **a native-free consensus/disagreement
signal is anti-informative about whether fusing helps**, which sits directly against the standing
project law that consensus is the only in-band discriminator, and against the 3D-Jury /
consensus-QA premise that pairwise similarity predicts quality.

---

# TASK B — SHAO & ZHU (2018), REINSTATED AND RE-SCOPED

## B.1 The source, obtained

RSC returns 403 and PubMed a cookie wall, as LIT reported. The **publisher-deposited abstract record
was obtained verbatim from the Europe PMC REST API**
(`/europepmc/webservices/rest/search?query=DOI:"10.1039/c7cp08010g"&resultType=core&format=json`).
**Tier [A+]:** the deposited record itself, not a search snippet and not a summarising model's
rendering of a page — a strictly better tier than anything LIT could reach. **Full text NOT
obtained** (`isOpenAccess: N`, `inEPMC: N`, subscription only), so anything below Methods level is
inference and is labelled.

> Shao Q, Zhu W. "Assessing AMBER force fields for protein folding in an implicit solvent."
> *Phys. Chem. Chem. Phys.* **20**(10):7206–7216 (2018). doi:10.1039/c7cp08010g. PMID 29480910.
> 39 citations.

**Verbatim:** *"we performed enhanced sampling MD simulations to assess the ability of six AMBER
force fields (FF99SBildn, FF99SBnmr, FF12SB, FF14ipq, **FF14SB**, and **FF14SBonlysc**) as coupled
with a recently improved pair-wise **GB-Neck2** model in modeling the folding of two helical and two
β-sheet peptides. Whilst most of the tested force fields can yield roughly similar features for
equilibrium conformational ensembles and detailed folding free-energy profiles for short α-helical
TC10b in an implicit solvent, the measured counterparts are significantly discrepant in the cases of
larger or β-structured peptides (HP35, 1E0Q, and GTT). Additionally, the calculated folding/unfolding
thermodynamic quantities can only partially match the experimental data. Although a combination of
the force fields and GB-Neck2 implicit model able to describe all aspects of the folding transitions
towards the native structures of all the considered peptides was not identified, we found that
FF14SBonlysc coupled with the GB-Neck2 model seems to be a reasonably balanced combination to predict
peptide folding preferences."*

**Systems, lengths pinned against the RCSB entry API (`deposited_polymer_monomer_count`):**

| system | residues | class | verdict |
|---|---|---|---|
| TC10b (Trp-cage variant, cf. 2JOF) | **20** | α | works |
| HP35 (villin headpiece, cf. 1YRF) | **35** | α (3-helix) | discrepant |
| **1E0Q** (N-terminal 17-mer of ubiquitin) | **17** | β-hairpin | **discrepant** |
| GTT (FiP35 variant, Pin1 WW domain) | ~35 | β | discrepant |

TC10b↔2JOF and GTT↔FiP35 identifications are **[I]**; the 1E0Q entry title is *"Mutant Peptide from
the first N-terminal 17 amino-acid of Ubiquitin"*, 17 residues, and is **[verified]**.

**LIT's condition is confirmed: this is our exact ff14SB/GB-Neck2 pair, and Sprint 15's stated
ground for the downgrade fails.**

## B.2 A correction to LIT, in Sprint 15's favour

LIT wrote: *"Sprint 15 §2.4 downgraded the threat ... on the stated grounds that it does not test
ff14SB."* **That conflates two sources.** Sprint 15 downgraded **Maffucci & Contini (2016)**, for
which the ground is true, and separately recorded the *PCCP* study as **UNVERIFIED**, flagged it as
*"the highest-value unfetched source in this document"*, and pre-registered the exact consequence:
*"Do it before any paper draft cites AMBER anti-ranking as novel."* **Sprint 15's judgement here was
correct and its risk was correctly localised.** What is reversed is the *combined* §2.4 verdict, not
that bullet. Recorded because the programme's standard is that corrections run in both directions.

LIT also over-reads the conclusion. *"Finds no combination describes folding to the native"* is
stronger than the abstract, which says no combination described **all aspects** of the folding
transitions of **all** the peptides — while reporting the α-helical 20-mer as working and endorsing
**FF14SBonlysc/GB-Neck2 as "a reasonably balanced combination"**. It is a mixed result, and a
referee citing it against us is citing a mixed result.

## B.3 Is the re-scoping defensible or self-serving? — **hostile assessment**

LIT proposes: their measure is MD folding **thermodynamics**, ours is certified-optimum conformer
**ranking**. **Two of the three natural supports for that re-scoping fail, and they must not be
used.**

**FAILS — "they tested longer peptides than ours."** `1E0Q` is a **17-residue β-hairpin**, one
residue above our 9–16 band, and it is named among the discrepant systems. The length defence is not
available for the β-sheet failure. (It *is* available for HP35 and GTT at ~35 residues, and the one
system in their set closest to our regime *by class* is the one that fails.)

**FAILS, BACKWARDS — "free energy is a weaker observable than potential-energy ranking."** A folding
free-energy profile is Boltzmann-weighted and carries the conformational entropy our single-point and
minimised-pool rankings discard. **Their observable is the stronger one.** If ff14SB/GB-Neck2 does
not place the native basin at the free-energy minimum, a potential-energy ranking of a pool drawn
from the same landscape has no reason to do better. **Their result is upstream of ours and predicts
it.** Any framing that implies we tested something harder is self-serving and must be struck.

**HOLDS — the observable is genuinely different, and it is the only defensible part.** Nothing in the
paper reports a **rank position of a native within a fixed candidate pool**; there is no enumerated
space, no certified global optimum, no per-target confidence interval, no filtering/ordering
decomposition, and no polarizable comparison. n = 4 peptides against our 126 with pinned folds and 19
with a certified optimum over 1.28 × 10⁷ exactly labelled structures.

> **Verdict: the re-scoping is DEFENSIBLE ONLY IN ITS NARROWEST FORM, and Sprint 15's framing was
> too generous to itself.** The honest sentence is: *"ff14SB with GB-Neck2 is already documented not
> to describe folding to the native for β-structured peptides down to 17 residues (Shao & Zhu 2018).
> We do not discover that AMBER mis-ranks peptide natives; we quantify it in a different observable
> — rank position within a fixed pool against a certified global optimum — with paired fold-aware
> intervals on 126 targets, and we show it survives polarizable physics."*

## B.4 Reconciliation with what this programme has already measured

Every one of our own AMBER results is **consistent with** Shao & Zhu, which is the strongest reason
to treat the re-scoping as honest rather than convenient — and also the reason the discovery claim
must go.

| our record | relation to Shao & Zhu |
|---|---|
| **native anti-ranking is distributed, not a single-term defect** — a property of the total energy | Their shape exactly: *no combination* fixes *all aspects*. Neither result localises the failure to a removable term |
| **AMOEBA does not fix it** — polarization nearly flips trpzip, chignolin fails on vdW | **Genuinely beyond their scope.** They test six fixed-charge AMBER variants; a polarizable extension of the ceiling is ours. **This is the strongest surviving piece** |
| **converged interaction-only AMBER flips the sign on the decoy bank but does not transfer to real pools** | Consistent: garbage rejection is not nativeness. Their peptides are physical throughout, so they never see the decoy-bank effect at all |
| **on matched pools with an r_g control AMBER is significantly *worse* than the control and puts the native at the 51st percentile**; all-atom reranking inside a real pool is worth **+0.004 Å** | This is the *ranking* observable they do not measure. It is the piece that survives — and it says something sharper than "the thermodynamics are discrepant": within a realistic pool the energy carries **no** ordering information at all |
| **Sprint 15 PHYS: the −0.022 Å accuracy claim is ABSENT on the frame-reproducible 69% of the instrument** (−0.009 [−0.024, +0.006], 38W/46L, ×3 frames) | **Reinforces the reinstatement.** Our own strongest AMBER *accuracy* claim already failed its own controls. Nothing is lost by conceding the anti-ranking priority, because the accuracy half was already weakened internally |
| **Sprint 15 PHYS: the VALIDITY claim is CONFIRMED and larger than recorded** — Ramachandran 0.466 → **0.874** (116W/2L), clashes 1.397 → **0.000** (63W/0L) | **Untouched by this literature.** Shao & Zhu is about where the free-energy minimum sits, not about whether a restrained minimisation repairs stereochemistry. **This is the programme's one unthreatened positive physics result and it should carry the AMBER pillar** |

## B.5 Where the corrected scoping was written

| file | edit |
|---|---|
| `s15/LITERATURE.md` §2.4 | full correction box at the head: verbatim abstract, the four systems with residue counts, the hostile assessment, the corrected position. Original heading and text struck through, preserved |
| `s15/LITERATURE.md` §2.4 bullets | the "does NOT close our claim 3" bullet struck through with the reversal; the "Still UNVERIFIED" bullet marked RESOLVED, **with credit to Sprint 15 for localising the risk correctly** |
| `s15/LITERATURE.md` table row **C4b** | rewritten from *"[S] snippet only, UNVERIFIED"* to the full **[A+]** entry with the re-scoping |
| `s15/LITERATURE.md` §5.4 | the "RESOLVED IN OUR FAVOUR" bullet struck through and replaced with "NOT RESOLVED IN OUR FAVOUR" |
| `s15/LITERATURE.md` §5.5 claim table, row 3 | narrowed: the force-field half of the claim is also taken |
| `s15/LEGACY_VS_AMBER.md` §3.1 | Shao & Zhu added as published context beside Roget et al., with the two failed defences named first |

---

## Reproducibility, compute, and rules

* **Bit-exact reproduction** of `s15/coherence.py`'s two fits: max |Δ| = 0.000e+00 on all four
  quantities, 126/126. The determinism check ran in a **fresh interpreter**, per BRIEF §3.7.
* **`stable_rng` only.** No `hash()`. Seeds: `("retract_law","foldboot17")`,
  `("retract","lemma2")`, and `coherence.py`'s own `stable_rng(pdb,"coherence")` for the refits.
* **Statistics.** TARGET as the unit throughout; n printed; median beside mean; W/L; per-fold detail;
  both an i.i.d.-target and a fold-clustered interval. The A.2 result is **W/L 126/0 with 5/5 folds
  the same sign**, so the BRIEF §3.10 concentration trap does not apply; no near-even W/L with a CI
  excluding zero appears in this document.
* **Small-n.** Nothing here rests on a small-n read. The n = 8 smoke of `retract_exact` agreed with
  n = 126 on the *exactness* question (3.4e−15 vs 2.7e−15) — as it must, since exactness is a
  machine-precision question and not an effect size — but its *effect sizes* differed materially
  (medoid-frame penalty 0.117 at n = 8 against 0.174 at n = 126, `kabsch_buys_on_s` 0.102 against
  0.368). Both are reported; only the n = 126 numbers are quoted.
* **Benchmark.** The sealed 60-target set was not read, probed, or referenced.
* **Leakage.** Natives are read only to audit a published number. No quantity here enters a
  predictor. Every native-reading arm is labelled ORACLE DIAGNOSTIC.
* **Compute.** One heavy process at a time. `retract_exact` is 322 s single-process with
  `OMP/MKL/OPENBLAS_NUM_THREADS=2`, checkpointed every 10 targets. `retract_law` and
  `retract_screen` are JSON reads. Load was sampled (94–100% CPU, 4.1–4.6 GB free) before each
  launch; nothing was launched in parallel with anything else of mine.

## Recommended ledger entries for the coordinator

1. **The fusion law is withdrawn as a contribution** (prior art, wrong mean, exact identity, refuted
   screen). `s16/BRIEF.md` §2 still asserts it and should be corrected by its owner.
2. **A transcription audit cannot see a formula error.** The Sprint 15 verification pass validated
   this number to the artefact and the artefact was wrong. Any future audit pass needs an arm that
   re-derives a published identity from its definition.
3. **The Shao & Zhu downgrade is reversed**; "AMBER anti-ranks peptide natives" moves from finding to
   confirmation, and the AMBER pillar's defensible claims are the **observable**, the **statistics**,
   the **AMOEBA extension**, and the **stereochemical repair** result.

## STATUS: COMPLETE.
