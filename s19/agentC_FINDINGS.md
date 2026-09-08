# SPRINT 19 — AGENT C (LEGACY / AMBER / PHYSICS) — FINDINGS

**Genuine Legacy (`core.energy`, eleven components, `DEFAULT_WEIGHTS`, never fitted) and genuine
AMBER ff14SB/GBn2 (`core.amber`, OpenMM, CPU platform, `Threads=1`, `DeterministicForces`). No
learned surrogate stands in for either anywhere in this workstream.**

Pre-registration: `s19/PREREG_C.md`, written before any Sprint-19 physics number existed, with my
honest prior attached to every hypothesis — including my prior that the headline one would only
half work.
Modules: `s19/agentC_lib.py` (the common frame, the Krogh–Vedelsby decomposition, the gates, the
concentration check, gates GC0/GC1), `s19/agentC_kv.py` (the pre-registered primary),
`s19/agentC_kv2.py` and `s19/agentC_kv3.py` (two **declared extensions**, each labelled with what
it was written after), `s19/agentC_pareto.py` (Q2, gate GC2), `s19/agentC_reject.py` (Q3),
`s19/agentC_report.py` (every table). Artefacts in `s19/results/`, each with an explicit
`complete` flag, a row count and a config hash.

Claim labels: **EXACT** (a theorem or identity — not a discovery) · **ORACLE** (needs the native) ·
**ESTABLISHED · SUPPORTED · PLAUSIBLE · OPEN · INCONCLUSIVE · REFUTED**.

---

# 0. VERDICT — leading with what most damages my own hypothesis

> ## 1. MY OWN PRE-REGISTERED MECHANISM IS REFUTED BY A CONTROL I BUILT TO TEST IT, AND THE REFUTATION IS CLEAN.
>
> The pre-registered hypothesis (H-C1) was that a physics score gate damages the averaged answer
> because it **collapses the ensemble's ambiguity `D`** — the survivors' errors stop cancelling.
> In the exact Krogh–Vedelsby accounting that is spectacularly true: Legacy's gate loses
> **−Δ(D²) = +0.692 Å² [+0.581, +0.843], 14W/112L**, and the diversity channel accounts for
> **148%** of the damage while the member-quality channel points the *other* way.
>
> **And it is not the cause.** `rand_lowD` — the lowest-diversity of 200 matched-random subsets of
> the same size, containing **no score, no physics and no sequence information** — reaches
> **D = 1.723, which is Legacy's D to three decimals**, and costs
> **+0.003 Å [−0.003, +0.013]** where Legacy costs **+0.063 Å [+0.010, +0.124]**.
> The pure diversity contrast `rand_lowD − rand_highD` (D 1.723 vs 2.084, a *wider* gap than
> Legacy vs random) is **+0.008 Å [−0.013, +0.028] — a null.**
>
> **An accounting identity absorbed the damage; it did not cause it.** I report this as the
> primary result of my own arm because it is the one that most damages it.

> ## 2. THE MECHANISM IS NOT DIVERSITY. IT IS *WHERE* THE GATE MOVES THE ANSWER, AND IT SPLITS EXACTLY IN TWO.
>
> In one fixed common frame, with `b_pool` the pool's mean error and `Δ_S = b_S − b_pool` the
> displacement a gate applies to the emitted mean, the readout obeys, **exactly**:
>
>     readout_S² = ‖b_pool‖² + 2‖b_pool‖·ALIGN_S + ‖Δ_S‖²,     ALIGN_S := ⟨b_pool, Δ_S⟩/‖b_pool‖
>
> so a gate's damage has **a MAGNITUDE term and an ALIGNMENT term and nothing else**. A
> matched-random gate has `E[ALIGN] = 0` *whatever its diversity*, and pays only the magnitude
> term — which is the recorded **+0.013 Å cost of halving**. Measured, n = 126:
>
> | | ALIGN (Å) | ‖Δ‖ (Å) | damage vs matched-random |
> |---|---|---|---|
> | the five **physics score** gates, pooled | **+0.0209 [+0.0004, +0.0402]** | 0.365–0.478 | +0.023 to +0.063 |
> | the four **score-free** gates, pooled | **−0.0048 [−0.0124, +0.0020]** | 0.245–0.625 | −0.005 to +0.043 |
> | score-free **minus** score | **−0.0256 [−0.0455, −0.0008]** | | |
>
> **A physics score gate systematically pushes the emitted mean ALONG the direction in which the
> pool is already wrong. Score-free gates of the same count, at every diversity from the lowest to
> the highest available, do not.** (`helix` is in neither pool: it *is* a score, just a
> zero-information one, and it has the largest ALIGN of all — see §0.3. The dichotomy that matters
> is **ordering by a score** versus **not ordering by a score**, not physics versus not-physics.) Across 12 arms mean ALIGN predicts mean damage at
> **r = +0.838**; *within* an arm, per target, ALIGN excess predicts the damage at
> **r = +0.951 (legacy), +0.967 (helix), +0.955 (leg_steric)**, n = 126 each.
>
> Both pooled statistics pass a drop-top-vs-uniform-effect concentration check (null percentile
> 0.495 and 0.496 — **not concentrated**), which the near-even W/L makes mandatory.

> ## 3. THE ZERO-INFORMATION CONTROL IS THE WORST GATE IN THE TABLE, AND THAT IS THE FINDING, NOT AN ASIDE.
>
> `helix` — keep the 38 candidates closest to a **constant ideal α-helix**, a gate with no
> sequence, no structure prediction and no physics in it — costs **+0.1253 Å [+0.0563, +0.1819]**
> against matched-random, **twice the Legacy gate's +0.0631** and the largest ALIGN in the table
> (+0.0476 Å).
>
> **Physics scores are not being punished for being physics.** They are being punished for being
> *any* rule that concentrates the survivor set toward a common structure, and a rule with zero
> information does it hardest. This kills the reading of Sprint 18's hardest negative as
> "the physics is anti-correlated with quality"; the correct reading is
> **"an ordering — any ordering — displaces the emitted mean, and orderings that correlate with
> the pool's error displace it the wrong way."**

> ## 4. THE DESIGN RULE THAT FOLLOWS WORKS, AND IT RECOVERS THE WHOLE SPRINT-18 LEGACY DEFICIT — AND IT IS STILL NOT A DEPLOYABLE GAIN.
>
> `legacy_clust` (partition the 75 into 38 Cα-RMSD clusters, keep each cluster's **best-Legacy**
> member — same score, same count, native-free) against the plain Legacy gate:
> **−0.0664 Å [−0.1392, −0.0087], 75W/51L.** `legacy_spread` gives **−0.0610 [−0.1304, −0.0012]**.
> Both have **negative** ALIGN (−0.016) where plain Legacy is +0.026.
>
> **H-C2's first half is SUPPORTED and its second half is not.** Against a matched-random gate
> `legacy_clust` is **−0.0033 Å [−0.0198, +0.0112]**, SE 0.0174, **per-comparison MDE 0.049 Å** —
> so this comparison detects 0.049 Å and finds nothing. It **excludes any effect larger than
> ~0.05 Å** and cannot speak below that, which is tighter than the blanket 0.084 Å I first quoted
> (§1.2). It is a null with a stated resolution, not a "match". Restructuring a physics gate so that it stops moving the answer the wrong way buys
> back everything the ordering cost **and then stops exactly at the random gate.** The value of
> the physics inside it remains zero.

> ## 5. Q2 DELIVERED, AND F-C2 FIRES: NO CONSTRAINED AMBER BUYS THE REPAIR FOR LESS DISPLACEMENT.
>
> Sixteen protocols, three orders of magnitude of restraint strength, three restrained atom sets,
> n = 126, `steps = 0` (the deployed protocol), the full pre-registered ladder. **The frontier is a
> single monotone trade between displacement and repair and the incumbent k = 30 sits on it.**
> Every arm whose validity matches k30 is within **±0.024 Å** of it on accuracy — three of those
> below the 0.0142 Å rotated-frame floor, hence NOT MEASURED. The three arms that beat k30 past the
> MDE do it by **not repairing**: `blend75` −0.119 Å with **+11.63 extra close contacts**,
> `blend50` −0.092 Å with **+7.49**, and `caonly_k300` −0.110 Å with **+0.424 cis peptide bonds**.
> Sprint 18's *"the repair tax is irreducible"* reproduces by a completely different route.
>
> **And the sharpest single item in my lane is a methodological one.** `caonly_k300` is
> **0.110 Å [0.092, 0.136] more accurate than the incumbent with a clash count statistically
> indistinguishable from it** (−0.008 [−0.050, +0.032]). **A scalar validity score built on clashes
> would have promoted it as the sprint's headline.** What kills it is `cis_frac` **+0.4242**, and
> the mechanism is confirmed by its dose-response — pinning Cα harder runs cis 0.010 → 0.051 →
> 0.101 at k = 30 → 100 → 300 while the clash count stays at zero. Pin only Cα and the peptide
> plane flips. **One omitted axis and the conclusion inverts.**

> ## 5b. AND I PUT A FALSE DEFECT INTO THE PROGRAMME'S LEDGER, THEN REFUTED IT MYSELF.
>
> I claimed `core.amber`'s minimisation does not terminate on some inputs (">40 min at k = 3 on
> 1A1P"). **Re-measured on a quiet box, uncapped, the same calls take 14.7 s, 15.2 s and 12.8 s and
> all converge.** What I measured was **CPU starvation** from six other lanes on a shared box, and I
> attributed a scheduling artefact to the minimiser. There is no pathology. Three pre-registration
> deviations that rested on it are **withdrawn** and the run above is the pre-registration as
> written. **And the gate I built to police the fix had a second flaw the first one hid**: GC3
> certified the bound bit-inert on five targets where the incumbent converges *before* the bound —
> the one call the bound actually bound moved the structure **0.16 Å** and was never checked. Both
> are in §4.8–4.9, with the general rule: **an inertness certificate must be evaluated on the calls
> where the intervention acted.**

> ## 5c. AND I APPLIED THE BRIEF'S MDE CONSTANT AS A THRESHOLD, WHICH IT IS NOT.
>
> `MDE = 2.8016·SE` is a property of a **comparison**, not of an instrument. Across this lane the
> per-comparison MDE runs **0.010 → 0.113 Å, an 11× spread**, so the blanket 0.084 Å was too strict
> in Q2 (it hid **four** genuinely measured arms, taking the count that beat k30 from three to
> **six** — all six still failing validity) and too loose in Q1 (it made a 0.049 Å-resolution null
> look uninformative). **No verdict reverses**; F-C1, F-C2 and F-C3 all fire on the corrected
> arithmetic. The exposed caveat is real and I state it rather than bury it: **most of this lane's
> positive effects sit at 0.7–1.3× their own MDE** — significant by CI, underpowered by design, so
> their **magnitudes are upper bounds**. The directions reproduce Sprint 18 to the third decimal.
> The one place I have power to spare is the **null** side that refutes my own hypothesis. §1.2.

> ## 6. THE OPERATOR LAW'S FAILURE IS NOW EXPLAINED, AND IT REPRODUCES TO THREE DECIMALS.
>
> Recomputed on my own arms, not imported: the miss of `d_out = 1.16·d_set_mean + 0.04·d_set_best`
> is **+0.0937 [+0.0629, +0.1317]** for Legacy (Sprint 18: +0.094 [+0.063, +0.132]) and
> **+0.0790 [+0.0601, +0.0983]** for `leg_torsion` (Sprint 18: +0.079 [+0.060, +0.098]), against
> **+0.0098 [+0.0071, +0.0132]** for a matched-random gate. **The law consumes only `⟨d⟩` and
> `min d`. It is blind to `ALIGN` and to `‖Δ‖`, which is where all the damage is** — and for a
> random gate `ALIGN` is zero in expectation, which is exactly why the law holds there and only
> there. The law is not wrong; it is **incomplete in a way that is invisible until you intervene
> selectively.**

---

## CLAIM TABLE

| claim | tier | evidence |
|---|---|---|
| `readout² = E_mem² − D²` in the operator's common frame | **EXACT** (gate GC0) | max residual **8.5e-14 Å²** over 10 arms × 126 targets; `Y.mean(0)` reproduces the deployed average to **0.00e+00 Å** |
| `readout² = ‖b_pool‖² + 2‖b_pool‖·ALIGN + ‖Δ‖²` | **EXACT** | algebra; both terms measured separately in `agentC_kv3.json` |
| a physics gate **collapses ensemble diversity** relative to matched-random | **ESTABLISHED** | Legacy −Δ(D²) **+0.692 [+0.581, +0.843], 14W/112L**; every score gate except `leg_contact` the same sign |
| …and **that is why** the answer degrades — **my pre-registered H-C1** | **REFUTED, by my own control** | score-free `rand_lowD` at **the same D (1.723)** costs +0.003 [−0.003, +0.013] vs Legacy's +0.063; `rand_lowD − rand_highD` = +0.008 [−0.013, +0.028] |
| the damage is carried by **ALIGN**, the displacement along the pool's own error direction | **SUPPORTED** | score gates +0.0209 [+0.0004, +0.0402] vs score-free −0.0048 [−0.0124, +0.0020]; difference −0.0256 [−0.0455, −0.0008]; cross-arm r = +0.838, within-arm per-target r = +0.951 |
| …and the ALIGN effect is **not carried by a few targets** | **ESTABLISHED** | drop-top-5 vs uniform-effect null: percentile **0.495**, verdict *not concentrated* |
| Sprint 18's "every physics filter loses to a random gate" | **ESTABLISHED, reproduced with a 10× tighter null** | legacy **+0.0631 [+0.0097, +0.1236]** (S18 +0.063 [+0.008, +0.125]); `leg_torsion` +0.0415 (S18 +0.041); `leg_contact` +0.0534 (S18 +0.053); `amber_sp` +0.0363 (S18 +0.036, CI now spans zero at 30 draws); halving alone +0.0132 (S18 +0.013) |
| a **zero-information** gate is the *worst* of all | **ESTABLISHED — and it reframes the Sprint-18 negative** | constant-α-helix gate **+0.1253 [+0.0563, +0.1819]**, 2× Legacy's cost, largest ALIGN (+0.0476) |
| `leg_contact` damages by the **same** mechanism as the others | **REFUTED — it is the one exception** | its −Δ(D²) is **+0.009 [−0.147, +0.127]** (a null) and **96%** of its damage is member quality, Δ⟨d_reb²⟩ +0.398 [+0.051, +0.719]. It is a genuine anti-ranker; the others are not |
| a **diversity-preserving** Legacy gate beats the plain Legacy gate | **SUPPORTED (H-C2, first half)** | `legacy_clust` −0.0664 [−0.1392, −0.0087] 75W/51L; `legacy_spread` −0.0610 [−0.1304, −0.0012] |
| …and beats a **matched-random** gate | **REFUTED at this comparison's resolution (H-C2, second half)** | −0.0033 [−0.0198, +0.0112], SE 0.0174, **own MDE 0.049 Å** — detects 0.049 Å, finds nothing; silent below that |
| the ledger's operator law survives a **selective** intervention | **REFUTED, reproduced** | miss +0.0937 (legacy) / +0.0790 (`leg_torsion`) vs +0.0098 (random) |
| my Legacy scores **are** the Sprint-18 artefact's | **EXACT** (gate GC1) | max \|Δ term\| **0.00e+00** and max \|Δ ORACLE d\| **0.00e+00 Å**, 8 targets × 75 candidates × 11 terms |
| `_run_ref` **is** `core.amber._run` and the patched restraint set takes | **EXACT** (gate GC2) | 0.00e+00 Å / 0.00e+00 kcal against `_run`; restraint sizes 42 / 14 / 42; the third leg against the deployed `refine_coords` also gave 0.00e+00 Å in an earlier logged run and was recorded as **SKIPPED** (memory ceiling) in the persisted artefact — see the note under §1 |
| the brief's **0.084 Å MDE** is a valid significance threshold | **REFUTED — MDE is per-comparison** (§1.2) | `MDE = 2.8016·SE`; across this lane it runs **0.010–0.113 Å**. Too strict in Q2 (hid 4 measured arms), too loose in Q1 (a 0.049 Å-resolution null read as uninformative). **No verdict reverses**, but most positive effects here sit at 0.7–1.3× their own MDE, so magnitudes are upper bounds |
| an inertness gate that **passes** has certified something | **REFUTED — a gate can pass VACUOUSLY** (§1.1, §4.9) | GC3 certified the iteration bound bit-inert at **0.00e+00 Å** at two different values, on five targets that converge *before* the bound — i.e. on calls the bound never touched. The one call it actually bound (`k30` on `8T61`) moved **0.1595 Å / 0.418 kcal-mol** and the gate had certified the opposite. **A certificate must be evaluated where the intervention acted**, and the number of firings reported |
| `core.amber`'s minimisation does not terminate on some inputs | **RETRACTED — my own claim, my own error** (§4.8) | re-measured uncapped on a quiet box: 1A1P k = 1/3/10 take **14.7 / 15.2 / 12.8 s** and converge; 1CEK ~6 s; 1D0W 7–10 s. Over ~2000 minimisations mean per-call wall **6.5 s**, p95 **13.6 s**. The original observation was **CPU starvation on a shared box**, misattributed to the minimiser |
| an inertness certificate computed where the intervention did not act is weak evidence | **REFUTED — it is NO evidence** (§4.9) | GC3 certified the bound at 0.00e+00 Å on five targets that converge *before* it; the one call it actually bound (`k30` on `8T61`) moved **0.1595 Å / 0.418 kcal-mol** and was never checked |
| a differently-constrained AMBER reaches a better Pareto point (Q2) | **REFUTED — F-C2 fires** (§4.4) | 16 protocols, n = 126: every arm whose validity matches k30 is within ±0.024 Å of it; the three that beat it past the MDE carry +11.63 / +7.49 extra clashes or +0.424 cis bonds |
| AMBER k = 30 is genuine stereochemical repair | **ESTABLISHED, sixth reproduction** (§4.3) | input clash 2.6 Å **1.113 → 0.000**, bond_strain **0.131 → 0.011**, min_heavy 2.752 → 2.886 |
| a **scalar** validity score would have reached the same verdict | **REFUTED — the sprint's cleanest vindication of the vector rule** (§4.5) | `caonly_k300` is **−0.110 Å [−0.136, −0.092]** on Cα with clash count indistinguishable from k30 (−0.008 [−0.050, +0.032]); only `cis_frac` **+0.4242** disqualifies it, dose-dependent 0.010 → 0.051 → 0.101 in the Cα pin |
| the AMBER convergence gate excludes the same three targets | **ESTABLISHED, 6th confirmation** (§4.1) | `1D6X 2NB7 7BX2`; and exclusions rise **3 → 24 → 64** at k = 100 → 300 → 1000, because a pinned backbone cannot relieve its own strain |

---

# 1. GATES — passed before any scientific number was quoted

`python -m s19.agentC_lib` → `s19/results/agentC_gates.json`;
`python -m s19.agentC_pareto --gate`.

| gate | what it checks | result |
|---|---|---|
| **GC0** | `common_frame_members` **is** the deployed averaging operator (`s15.phys_repl.averaged_backbone_from`); its mean **is** the metric's argument; the KV identity holds; `FRAME² ≥ 0` | **PASS** — mean delta **0.00e+00 Å**, readout-vs-metric 8.1e-15 Å, KV residual **3.4e-14 Å²**, min FRAME² **+5.2e-04 Å²** |
| **GC1** | my Legacy components and ORACLE labels **are** `s18/results/down.json`'s, candidate by candidate | **PASS** — **0.00e+00** on both |
| **GC2** | `_run_ref(H, pos, pos, k)` **is** `core.amber._run(H, heavy, k)` **is** the deployed `refine_coords`; the patched restraint set really changes `_restraint_idx` and the module constant is restored | **PASS** — 0.00e+00 Å and 0.00e+00 kcal/mol on both comparisons; 42 → 14 restrained atoms; `RESTRAINED_BACKBONE` restored |
| **GC3** | the iteration bound is **inert**: capped k = 30 must reproduce uncapped k = 30 bit-for-bit | **PASSED at both values tried, AND THE CERTIFICATE WAS WORTHLESS BOTH TIMES.** 0.00e+00 Å over 5 targets at `steps = 10000` and again at `steps = 2000` — but those five targets converge **before** the bound, so the gate certified calls the bound never touched. The one call it did bound (`k30` on `8T61`) moved **0.1595 Å** and was never checked (§4.9). **The bound has since been withdrawn entirely** (§4.8); the deployed Q2 run uses `steps = 0` and needs no such gate |

**A discipline note on GC2's third leg.** GC2 compares `_run_ref` against `core.amber._run` **and**
against the deployed `refine_coords`. The third comparison goes through `core.amber.builder_for`,
whose `memory_guard()` refuses to build a Context above 92% physical memory, and on the
`agentC_gates_amber.json` run it raised `MemoryError: physical memory at 93%` and was recorded as
**SKIPPED** (`deployed_leg_run: false`) rather than silently passed. That leg **did** pass at
**0.00e+00 Å** in an earlier logged run of the same code path, before the iteration bound was added.
Both facts are stated; neither is glossed.

### 1.1 A GATE CAN PASS VACUOUSLY — the rule this lane produced, and it applies to my own gate

GC3 certified an iteration bound as **bit-inert at 0.00e+00 Å**, twice, at two different values —
and both certificates were worthless:

> **the bound was inert on the five targets tested because on those targets it was never reached.**

A "this intervention changes nothing" certificate means nothing unless the intervention was
**exercised**. GC3's five targets are ones where the incumbent k = 30 converges in 3–9 s, well
before any bound; so the gate certified calls the bound never touched, and said nothing about the
one call it did. Checked afterwards, directly, that call disagrees:

    8T61  k30  bounded(2000) vs unbounded:  max delta Ca = 1.595e-01 A,  delta E = 0.418 kcal/mol

**The bound bound, it moved the structure by 0.16 Å, and the gate had certified the opposite.**
An arm that hits a bound is a **different operator** — "N iterations of restrained relaxation", not
"restrained relaxation".

**The general form, for any future inertness certificate in this programme:** evaluate it **on the
calls where the intervention actually acted**, and report how many times the bound, guard, threshold
or fallback **fired**. A certificate computed only on the cases the intervention did not touch is
not weak evidence — it is *no* evidence, and it reads as a pass.

The bound has since been **withdrawn entirely** (§4.8): the defect that motivated it does not exist,
and the deployed Q2 run uses `steps = 0` with no bound anywhere.

### 1.2 THE 0.084 Å MDE CONSTANT IS NOT A THRESHOLD, AND I USED IT AS ONE

`s19/BRIEF.md` §1 states a single instrument-wide MDE of **0.084 Å**, and I applied it as a blanket
significance threshold in both H-C2 and F-C2. **That is wrong.** MDE is a property of *a comparison*,
not of an instrument: `MDE = 2.8016 · SE`, and SE depends on the variance of the paired difference,
which varies by an order of magnitude across the contrasts in this lane. Corrected, with SE reported
beside every decisive mean as it should have been from the start:

| comparison | n | mean (Å) | SE | own MDE | \|mean\|/MDE |
|---|---|---|---|---|---|
| **Q1** legacy gate vs matched-random | 126 | +0.0631 | 0.0257 | 0.072 | 0.88 |
| **Q1** zero-information helix vs matched-random | 126 | +0.1253 | 0.0402 | 0.113 | **1.11** |
| **Q1** cost of halving alone | 126 | +0.0132 | 0.0037 | 0.010 | **1.27** |
| **Q1** `legacy_clust` vs plain legacy (H-C2 ①) | 126 | −0.0664 | 0.0335 | 0.094 | 0.71 |
| **Q1** `legacy_clust` vs matched-random (H-C2 ②) | 126 | −0.0033 | 0.0174 | 0.049 | 0.07 |
| **Q1** `rand_lowD` vs matched-random (score-free, matched D) | 126 | +0.0029 | 0.0132 | 0.037 | 0.08 |
| **Q1** `rand_lowD` − `rand_highD` (pure diversity contrast) | 126 | +0.0075 | 0.0192 | 0.054 | 0.14 |
| **Q1** ALIGN excess, 5 physics score gates | 126 | +0.0209 | 0.0070 | 0.020 | **1.06** |
| **Q1** ALIGN excess, 4 score-free gates | 126 | −0.0048 | 0.0074 | 0.021 | 0.23 |
| **Q1** score-free minus score (ALIGN) | 126 | −0.0256 | 0.0104 | 0.029 | 0.88 |
| **Q3** `steric@10` vs matched-random | 126 | +0.0175 | 0.0067 | 0.019 | 0.93 |
| **Q3** `steric@10` vs no rejection | 126 | +0.0204 | 0.0074 | 0.021 | 0.98 |
| **Q2** `caonly_k300` vs k30 | 123 | −0.1100 | 0.0109 | 0.030 | **3.62** |
| **Q2** `k100` vs k30 | 123 | −0.0386 | 0.0057 | 0.016 | **2.42** |
| **Q2** `bbo_k30` vs k30 | 123 | +0.0106 | 0.0067 | 0.019 | 0.57 |

**The per-comparison MDE in this lane runs from 0.010 Å to 0.113 Å — an 11× spread.** The blanket
constant was too strict in Q2 (it hid four genuinely measured arms) and too loose in Q1 (it made a
tight null look uninformative).

**What changes, stated plainly:**

* **F-C2 is strengthened.** Under the constant, three arms beat k30 past "the MDE"; under their own
  MDEs, **six** do — and all six still fail the validity axis with CIs excluding zero. §4.4.
* **H-C2's second half is tightened, not weakened.** `legacy_clust` vs matched-random resolves
  **0.049 Å** and returns −0.003: it excludes any gain above ~0.05 Å rather than the ~0.084 Å I
  claimed.
* **What does NOT change: no verdict in this document reverses.** F-C1, F-C2 and F-C3 all fire on
  the corrected arithmetic.

**And the caveat this exposes, which I would rather state than have read off the table.** Most of
this lane's *positive* effects sit at **0.7–1.3 × their own MDE** — significant by a fold-clustered
CI, underpowered by design. That is precisely the regime in which observed effect sizes are inflated
(the winner's curse). **The directions are safe** — they reproduce Sprint 18 to the third decimal,
they are consistent across five scores and twelve gate designs, and Q3's dose-response is monotone
across four doses. **The magnitudes should be read as upper bounds.** The one Q1 statistic that is
well powered on the *null* side is the pair that carries the refutation of my own hypothesis
(`rand_lowD` vs random at 0.08 of its MDE, the pure diversity contrast at 0.14), which is the right
way round: **the claim I killed was killed with power to spare.**

---

### 1.3 A defect GC0 caught in its own first draft, reported rather than repaired silently

My first (KV3) used the operator law's `d` — the free-superposition RMSD of the **raw window CA
trace** — as the member error, and `FRAME² = E_mem² − ⟨d²⟩` came out **negative, −0.059 Å²**, which
is impossible for a minimum-over-rigid-motions gap. **Cause**: the operator averages the
*ideal-geometry rebuild* of each window's (φ, ψ), not the window. The two differ by the project's
0.347 Å rebuild floor. `agentC_lib` now carries both quantities under different names and never one
in the other's place. This is exactly the error class the brief names (a prior sprint broke the same
identity by 18%), and the gate caught it before it reached a result.

---

# 2. Q1 — WHY SCORE-ORDERING DAMAGES AN AVERAGED SET

Artefacts: `s19/results/agentC_kv.json` (pre-registered primary, complete, 126 rows),
`agentC_kv2.json`, `agentC_kv3.json` (declared extensions, complete, 126 rows each).
Tables: `s19/results/agentC_kv_report.txt`, `agentC_kv2_report.txt`, `agentC_kv3_report.txt`.

## 2.1 The instrument, and what in it is a theorem

Every arm keeps **m = 38 of the same K = 75 shipped windows** and hands them to the deployed
all-atom coordinate-average operator. In that operator's own common frame the parallel-axis
theorem gives, exactly,

    readout² = E_mem² − D²          and       readout² = ⟨d_reb²⟩ + FRAME² − D²

verified to **8.5e-14 Å²** on every arm and target. **Nothing in that sentence is a discovery.**
The scientific content is entirely in *how a gate moves the right-hand terms*.

## 2.2 The levels (n = 126)

| arm | readout | E_mem | D | FRAME | err_cos | shared |
|---|---|---|---|---|---|---|
| `none` (m = 75) | 3.050 | 3.720 | 1.930 | 0.725 | 0.684 | 0.669 |
| `rand` (30 draws) | 3.063 | 3.720 | 1.913 | 0.731 | 0.685 | 0.675 |
| `legacy` | **3.126** | 3.678 | **1.723** | 0.690 | **0.732** | **0.724** |
| `leg_torsion` | 3.105 | 3.670 | 1.730 | 0.652 | 0.728 | 0.718 |
| `leg_contact` | 3.116 | 3.759 | 1.917 | 0.735 | 0.691 | 0.683 |
| `leg_steric` | 3.086 | 3.687 | 1.819 | 0.668 | 0.706 | 0.696 |
| `amber_sp` | 3.099 | 3.711 | 1.845 | 0.685 | 0.706 | 0.695 |
| `helix` **(zero-information)** | **3.188** | 3.658 | **1.573** | 0.649 | **0.769** | **0.766** |
| `legacy_spread` | 3.065 | 3.774 | 2.018 | 0.783 | 0.662 | 0.657 |
| `legacy_clust` | 3.060 | 3.765 | 2.009 | 0.772 | 0.663 | 0.658 |

`err_cos` is the mean pairwise cosine of the members' **error vectors** and `shared` is
`readout²/E_mem²`. Both say the same thing the diversity column says: a score gate leaves survivors
whose errors are more nearly parallel. **It is a real, large, reproducible correlation — and §2.4
shows it is not the cause.**

## 2.3 The pre-registered primary, reported exactly as specified

`PREREG_C.md` §3.4 required three conditions on the **Legacy** arm.

| condition | result | verdict |
|---|---|---|
| (1) `mean Δ(readout²) > 0`, fold CI excluding zero | **+0.4671 [−0.0308, +1.1057]**, med +0.1026, 53W/73L | **fails as literally written** |
| (2) `mean(−ΔD²) > 0`, fold CI excluding zero | **+0.6922 [+0.5808, +0.8426]**, med +0.4646, **14W/112L** | **holds, decisively** |
| (3) `S_D ≥ 0.50` | **S_D = 1.48** | **holds** |

**Condition (1) is a defect in my own pre-registration, and I am keeping it.** The damage is
unambiguous on the **metric** scale — `legacy` vs `rand` is **+0.0631 Å [+0.0097, +0.1236]**, which
reproduces Sprint 18's +0.063 [+0.008, +0.125] with a 10× tighter null — but I specified the
primary on the **squared** scale so that the channel shares would be exact, and the squared scale is
heavy-tailed enough that its CI spans zero at n = 126. Both numbers are reported side by side and
neither is substituted for the other.

### The full channel table (Å², gate minus matched-random, exact and additive)

| arm | Δ(readout²) | Δ(E_mem²) | −Δ(D²) | Δ⟨d_reb²⟩ | Δ(FRAME²) | S_D |
|---|---|---|---|---|---|---|
| `legacy` | +0.467 | −0.225 | **+0.692** | −0.147 | −0.078 | 1.48 |
| `leg_torsion` | +0.430 | −0.161 | **+0.591** | −0.012 | −0.149 | 1.38 |
| `leg_contact` | +0.414 | **+0.405** | +0.009 | **+0.398** | +0.007 | **0.02** |
| `leg_steric` | +0.113 | −0.311 | +0.424 | −0.159 | −0.152 | 3.74 |
| `amber_sp` | +0.265 | −0.012 | +0.278 | +0.099 | −0.111 | 1.05 |
| `helix` | **+0.956** | −0.185 | **+1.141** | −0.038 | −0.147 | 1.19 |
| `legacy_spread` | −0.118 | +0.324 | −0.441 | +0.212 | +0.112 | — |
| `legacy_clust` | −0.155 | +0.237 | −0.393 | +0.153 | +0.085 | — |

Two things to read off it. **First, `leg_contact` is the exception**: it is the only gate that does
not collapse diversity (−Δ(D²) = +0.009 [−0.147, +0.127], a null) and the only one whose damage is
**member quality** (Δ⟨d_reb²⟩ = +0.398 [+0.051, +0.719], CI excluding zero). That is consistent with
Sprint 17/18's measurement that `leg_contact` is a genuine *anti*-ranker (ρ_global −0.176); the
other four are not anti-rankers at all — **their survivors are measurably BETTER
(Δ(E_mem²) < 0) and the answer is still worse.** Second, on the four non-`leg_contact` gates the
member-quality and FRAME channels both point the *right* way and are overwhelmed by the diversity
term.

## 2.4 THE DECLARED EXTENSION THAT KILLS MY OWN MECHANISM

`agentC_kv2.py`, written after §2.3 was read and labelled as such. It builds gates with
**no score, no physics and no sequence information**, and manipulates diversity directly:
`rand_lowD` / `rand_highD` are the lowest- and highest-`D` of 200 matched-random subsets of the
same size (`D` is native-free, so this is a legitimate native-free construction); `fps` is pure
greedy farthest-point; `medoid` is the m candidates nearest the pool medoid.

| arm | D | readout | vs matched-random |
|---|---|---|---|
| `legacy` | **1.723** | 3.126 | **+0.0631 [+0.0097, +0.1236]** |
| `rand_lowD` **(score-free, same D)** | **1.723** | 3.066 | **+0.0029 [−0.0034, +0.0128]** |
| `rand_highD` (score-free) | 2.084 | 3.058 | −0.0046 [−0.0295, +0.0177] |
| `medoid` (score-free, extreme) | 1.419 | 3.106 | +0.0430 [+0.0131, +0.0762] |
| `fps` (score-free, max diversity) | 2.019 | 3.063 | −0.0003 [−0.0130, +0.0164] |
| `helix` (zero-information score) | 1.573 | 3.188 | +0.1253 [+0.0563, +0.1819] |

> **At the same diversity, the score gate costs 20× what the score-free gate costs.**
> `rand_lowD − rand_highD` — a *wider* diversity contrast than `legacy` vs `rand` — is
> **+0.0075 [−0.0128, +0.0279]**, a null.
>
> Pooled over 13 arms × 126 targets, `corr(Δreadout², −ΔD²) = +0.174` against
> `corr(Δreadout², ΔE_mem²) = +0.862`.

**H-C1 is REFUTED as a causal mechanism.** Diversity collapse is a real and large *companion* of
the damage — it is forced to absorb it by an identity — and is not what produces it.

## 2.5 THE SECOND DECLARED EXTENSION, AND THE MECHANISM THAT SURVIVES

`agentC_kv3.py`, written after §2.4. In one **fixed** common frame (the full pool's), with
`e_i = Y_i − T`, `b_pool = mean_K e_i` and `Δ_S = b_S − b_pool`:

    readout_S² = ‖b_pool‖² + 2‖b_pool‖·ALIGN_S + ‖Δ_S‖²          (EXACT)
    ALIGN_S = ⟨b_pool, Δ_S⟩ / ‖b_pool‖

`‖b_pool‖ = 3.050 Å` (= the ungated readout). `Δ_S` is **native-free**; `ALIGN` is ORACLE. The
matched-random ALIGN null, measured over 200 draws per target, is **+0.0008 Å, sd 0.0753**.

| arm | ALIGN (Å) | ‖Δ‖ (Å) | ALIGN term (Å²) | MAG term (Å²) | damage vs rand (Å) |
|---|---|---|---|---|---|
| `helix` | **+0.0476** | 0.694 | +0.305 | +0.624 | +0.1253 |
| `legacy` | +0.0256 | 0.478 | +0.187 | +0.263 | +0.0631 |
| `leg_contact` | +0.0169 | 0.460 | +0.183 | +0.225 | +0.0534 |
| `leg_torsion` | +0.0266 | 0.374 | +0.314 | +0.122 | +0.0415 |
| `medoid` (score-free) | **−0.0065** | **0.625** | −0.186* | **+0.479** | +0.0430 |
| `amber_sp` | +0.0258 | 0.367 | +0.188 | +0.127 | +0.0363 |
| `leg_steric` | +0.0136 | 0.365 | +0.065 | +0.123 | +0.0227 |
| `rand_lowD` (score-free) | +0.0052 | 0.302 | +0.134 | +0.058 | +0.0029 |
| `rand_highD` (score-free) | +0.0007 | 0.245 | −0.118 | +0.018 | −0.0046 |
| `fps` (score-free) | −0.0150 | 0.332 | −0.220 | +0.098 | −0.0003 |
| `legacy_spread` | −0.0157 | 0.329 | −0.230 | +0.096 | +0.0021 |
| `legacy_clust` | −0.0160 | 0.323 | −0.239 | +0.090 | −0.0033 |
| `rand` (the null) | +0.0008 | 0.214 | — | — | — |

\* the ALIGN and MAG columns are **excesses over matched-random**, so they sum to the total excess
in Å²; `medoid` is the instructive row — a score-free gate with the **second-largest displacement
in the table** and a **negative** ALIGN, whose damage is therefore **entirely** the magnitude term.

**The pooled test, which is the decisive one:**

    ALIGN excess, five PHYSICS SCORE gates    +0.0209 [+0.0004, +0.0402]  med +0.0074  54W/72L
    ALIGN excess, four SCORE-FREE gates       -0.0048 [-0.0124, +0.0020]  med +0.0009  63W/63L
    score-free MINUS score                    -0.0256 [-0.0455, -0.0008]  med -0.0004  63W/63L

    concentration check (drop-top-5 vs a UNIFORM-EFFECT null, mandatory beside a near-even W/L):
      ALIGN excess, score gates       mean +0.0209  med +0.0074  sd 0.0790  drop-top5 +0.0112
                                      null percentile 0.495  ->  NOT CONCENTRATED
      readout damage, legacy vs rand  mean +0.0631  med +0.0258  sd 0.2882  drop-top5 +0.0298
                                      null percentile 0.496  ->  NOT CONCENTRATED

**How strong is this, honestly.** The pooled score-gate ALIGN interval only just excludes zero
(+0.0004), and the per-arm intervals all span zero at n = 126 (`legacy` +0.0248 [−0.0157, +0.0738]).
What carries the claim is not any one interval, it is (a) the pooled score-vs-score-free contrast,
(b) the cross-arm ordering (r = +0.838 over 12 arms), and (c) the per-target correlation between
ALIGN excess and the realised damage: **+0.951, +0.967, +0.955** on `legacy`, `helix`, `leg_steric`
at n = 126. Label: **SUPPORTED**, not ESTABLISHED.

## 2.6 The design rule, and what it is worth

**The rule (SUPPORTED):** *a gate's damage is set by how far it moves the emitted mean and in which
direction — not by the set mean, not by the set best, and not by the diversity.* `‖Δ‖` is
**native-free and directly measurable**, so "how far did my gate move the answer" is a deployable
diagnostic; only the sign of the move (`ALIGN`) needs the native.

**What the rule buys (H-C2):** restructuring the Legacy gate so that it stops moving the answer the
wrong way recovers **the whole** Sprint-18 deficit —

| | vs plain `legacy` | vs matched-random | vs no gate |
|---|---|---|---|
| `legacy_clust` | **−0.0664 [−0.1392, −0.0087]** 75W/51L | −0.0033 [−0.0198, +0.0112] | +0.0100 [−0.0060, +0.0242] |
| `legacy_spread` | **−0.0610 [−0.1304, −0.0012]** 72W/54L | +0.0021 [−0.0089, +0.0124] | +0.0153 [+0.0050, +0.0250] |

— **and then stops at random.** Against matched-random the CI is far too wide to exclude the
effect of interest. Measured properly (§1.2): SE **0.0174**, **own MDE 0.049 Å** — the comparison
resolves 0.049 Å and returns −0.003, so it **excludes any gain larger than ~0.05 Å** and is silent
below that. Not a "match"; a null with a stated resolution.
The mechanism is confirmed; the deployable gain is zero. **Legacy's information content is still
worth nothing.**

## 2.7 The operator law, recomputed rather than imported

| arm | law predicts | observed | **miss [fold CI]** | Sprint 18 |
|---|---|---|---|---|
| `legacy` | −0.0173 | +0.0763 | **+0.0937 [+0.0629, +0.1317]** | +0.094 [+0.063, +0.132] |
| `leg_torsion` | −0.0243 | +0.0547 | **+0.0790 [+0.0601, +0.0983]** | +0.079 [+0.060, +0.098] |
| `helix` | +0.0031 | +0.1385 | **+0.1354 [+0.1048, +0.1643]** | — |
| `matched random` | +0.0034 | +0.0132 | **+0.0098 [+0.0071, +0.0132]** | +0.006 [−0.003, +0.014] |
| `legacy_clust` | +0.0504 | +0.0100 | **−0.0404 [−0.0599, −0.0231]** | — |

The law is blind to both terms of §2.5. For a matched-random gate the ALIGN term is zero in
expectation and the magnitude term is small, so the law works; for anything that orders the set it
does not, **and it can miss in either direction** — `legacy_clust` beats the law by −0.040.

---

## 2.8 The half of the mechanism that needs no native, and is therefore usable

`ALIGN` is ORACLE. **`‖Δ‖` is not** — it is just how far the gate moved the emitted mean away from
the *ungated* mean, and it can be computed for any candidate gate at design time with no native
anywhere.

    over 1512 (arm, target) pairs, 12 arms:
        corr(damage, |Delta|^2 excess)   = +0.169     NATIVE-FREE
        corr(damage, ALIGN term)         = +0.830     ORACLE
        corr(damage, the two summed)     = +0.841     (exact up to the square root)

    ranking the twelve GATE DESIGNS by their mean |Delta|:
        Spearman(mean |Delta|, mean damage) = +0.930          NATIVE-FREE

**Per target the native-free half is weak; per gate design it is strong.** So the usable rule is a
design-time one:

> **Rank candidate gates by how far they displace the emitted average from the ungated average.
> The one that moves it least is very likely the one that damages it least — ρ = 0.93 over twelve
> designs, with no native information.** Everything beyond that ranking (the *direction* of the
> move, which is what separates a physics gate from a score-free gate at the same displacement)
> requires the native and is therefore a diagnostic, not a selector.

This is stated as a **design rule with a measured native-free proxy**, not as an accuracy gain: the
best gate the rule selects is a random gate, and the incumbent already does not gate.

---

# 3. Q3 — CAN LEGACY REJECT A CLASS OF IMPOSSIBLE CANDIDATE? **NO. THE ROLE CLOSES.**

Artefact: `s19/results/agentC_reject.json` (complete, 126 rows).
Table: `s19/results/agentC_reject_report.txt`.

## 3.1 The class is real, and it is not a strawman

Before asking whether Legacy can reject impossible candidates, I measured how many there are. On the
ideal-geometry rebuild of the shipped top-75, per target:

    candidates with a heavy-atom contact < 2.0 A   mean 2.66 per 75,  median 1,  max 18,
                                                   57/126 targets have none
    candidates with a heavy-atom contact < 2.6 A   mean 20.1 per 75,  median 19,  max 55,
                                                   only 1/126 targets has none
    the worst contact in a pool                    mean 1.791 A,  minimum 0.104 A

**0.104 A between two heavy atoms is not a bad candidate, it is an impossible one**, and pools
contain them. So the question is well posed.

## 3.2 The pre-registered primary, r = 10

The rule declared in `PREREG_C.md` §5.3: the detector is a live role **iff it beats BOTH the
matched-random rejection AND the diversity-preserving rejection**, at r = 10, CI excluding zero.

| detector | vs `rand@10` (30 draws) | vs `divkeep@10` | vs no rejection |
|---|---|---|---|
| Legacy **`steric`** | **+0.0175 [+0.0048, +0.0348]** | **+0.0170 [+0.0091, +0.0286]** | **+0.0204 [+0.0081, +0.0376]** |
| `minheavy` (force-field-free) | **+0.0139 [+0.0025, +0.0238]** | +0.0133 [−0.0015, +0.0294] | **+0.0168 [+0.0044, +0.0278]** |
| Legacy total | +0.0182 [−0.0036, +0.0390] | +0.0176 [−0.0036, +0.0417] | +0.0211 [−0.0033, +0.0426] |
| genuine AMBER single point | +0.0109 [−0.0049, +0.0285] | +0.0103 [−0.0112, +0.0291] | +0.0138 [−0.0022, +0.0331] |

> **FALSIFIER F-C3 FIRES, and it fires with a sign.** `steric` does not merely fail to beat its
> controls — it is **significantly worse than both**, and significantly worse than not rejecting
> anything at all. An independent, force-field-free statement of the same physical claim
> (`minheavy`) fails the same way. **The role is CLOSED.**

## 3.3 The dose-response, which is the cleanest form of the Sprint-18 result

Mean Cα-RMSD of the emitted coordinate average, n = 126, as a function of how many of the 75 are
rejected:

| rule | r = 0 | r = 2 | r = 5 | r = 10 | r = 19 |
|---|---|---|---|---|---

# 4. Q2 — AMBER AS A CONSTRAINED TERMINAL OPERATOR: THE FRONTIER HAS NO BETTER POINT

Artefact: `s19/results/agentC_pareto.json` — **complete, 126 rows, `steps = 0` (the deployed
unbounded protocol), the FULL pre-registered ladder, 16 arms, ~2000 genuine ff14SB/GBn2
minimisations.** Table: `s19/results/agentC_pareto_report.txt`. Gates GC2, GC3 in
`agentC_gates_amber.json`.

> ## F-C2 FIRES. No constrained-AMBER protocol clears the pre-registered bar. The three arms that beat the incumbent on Cα by more than the MDE all do it by **not repairing**.

## 4.1 The convergence gate, declared before use — and it is itself a result

`PREREG_C.md` §2: converged iff the final energy, restraint off, is finite and ≤ 1000 kcal/mol.

| arm | excluded | which |
|---|---|---|
| `k1 k3 k10 k30 k100`, `caonly_k30/100/300`, `blend*` | **3** | **1D6X 2NB7 7BX2** |
| `stage_300_30`, `bbo_k30` | 2 | 1D6X 7BX2 |
| `k300` | 24 | |
| `pullback_30_1000` | 63 | |
| `k1000` | **64** | |

**The same three targets, for the sixth sprint running.** And the exclusion count is a physical
measurement, not bookkeeping: **the harder the restraint, the more targets end above 1000 kcal/mol
— 3/126 at k ≤ 100, 24 at k = 300, 64 at k = 1000.** A pinned backbone cannot relieve its own
strain. Hard restraints do not merely fail to help; **they fail to converge**, and any future arm
that "preserves the Cα trace" by pinning it hard inherits that.

This is why the Pareto test is **pairwise** (each arm against k30 on targets where both converge,
n = 102–123). An all-arm gate costs **64 of 126 targets to one divergent arm** and shifts the input
mean from **3.050 → 2.435 Å** — it silently selects the easy half of the instrument. Both are
printed; only the pairwise one is read.

## 4.2 AXIS 1 — accuracy, against the incumbent k = 30 (negative = better than k30)

| arm | ΔCα vs k30 [fold CI] | W/L | n |
|---|---|---|---|
| `k1` | **+0.0903 [+0.0574, +0.1293]** | 36/87 | 123 |
| `k3` | +0.0540 [+0.0325, +0.0734] | 34/89 | 123 |
| `k10` | +0.0238 [+0.0135, +0.0372] | 35/88 | 123 |
| `bbo_k30` | +0.0106 [−0.0024, +0.0286] | 56/67 | 123 |
| `stage_300_30` | +0.0078 [+0.0011, +0.0133] | 45/78 | 123 |
| `caonly_k30` | −0.0054 [−0.0107, +0.0008] | 63/60 | 123 |
| `k100` | −0.0386 [−0.0448, −0.0344] | 99/24 | 123 |
| `k1000` | −0.0493 [−0.0682, −0.0368] | 52/10 | 62 |
| `pullback_30_1000` | −0.0506 [−0.0744, −0.0378] | 53/10 | 63 |
| `blend25` | −0.0523 [−0.0625, −0.0452] | 112/11 | 123 |
| `k300` | −0.0737 [−0.0913, −0.0534] | 88/14 | 102 |
| `caonly_k100` | −0.0776 [−0.0925, −0.0666] | 93/30 | 123 |
| `blend50` | **−0.0922 [−0.1115, −0.0793]** | 107/16 | 123 |
| **`caonly_k300`** | **−0.1100 [−0.1355, −0.0924]** | 95/28 | 123 |
| `blend75` | **−0.1193 [−0.1458, −0.1018]** | 103/20 | 123 |

**The ladder is monotone in one thing only: how little the structure is allowed to move.** Mean Cα
displacement runs 0.497 → 0.420 → 0.352 → 0.302 → 0.237 → 0.169 → 0.098 Å from k = 1 to k = 1000,
and the accuracy penalty tracks it exactly. The input is the best Cα trace in the table and every
relaxation degrades it.

**The soft rungs are the interesting new data** (they are the two I wrongly dropped, §4.8):
`k1` is **+0.090 Å worse** than the incumbent and has **more Ramachandran outliers**
(+0.0340 [+0.0170, +0.0517]). Softening the restraint below k = 30 costs accuracy *and* buys no
validity. **There is no better Pareto point hiding at the soft end** — which means my unjustified
deviation cost the sprint nothing scientifically, though it was still unjustified.

## 4.3 AXIS 2 — physical validity, the same structures, never fused into a scalar

| | clash 2 Å | clash 2.6 Å | min_heavy | bond_strain | angle_strain | rama_fav | rama_out | cis_frac |
|---|---|---|---|---|---|---|---|---|
| **INPUT** | 0.177 | 1.113 | 2.752 | 0.131 | 0.059 | 0.974 | 0.011 | 0.000 |
| `k1` | 0.000 | 0.000 | 2.868 | 0.012 | 0.022 | 0.958 | 0.032 | 0.006 |
| `k10` | 0.000 | 0.000 | 2.882 | 0.010 | 0.019 | 0.982 | 0.008 | 0.002 |
| **`k30` (incumbent)** | 0.000 | **0.000** | 2.886 | **0.011** | 0.018 | 0.975 | 0.009 | 0.004 |
| `k100` | 0.000 | 0.016 | 2.891 | 0.018 | 0.020 | 0.977 | 0.014 | 0.000 |
| `k300` | 0.000 | 0.113 | 2.881 | 0.033 | 0.023 | 0.971 | 0.018 | 0.001 |
| `k1000` | 0.000 | 0.242 | 2.861 | 0.058 | 0.025 | 0.967 | 0.020 | 0.000 |
| `bbo_k30` | 0.000 | 0.000 | 2.874 | 0.010 | 0.018 | **0.990** | **0.004** | **0.000** |
| `caonly_k30` | 0.000 | 0.000 | 2.890 | 0.012 | 0.021 | 0.951 | 0.035 | 0.010 |
| `caonly_k100` | 0.000 | 0.000 | 2.890 | 0.021 | 0.024 | 0.927 | 0.058 | **0.051** |
| `caonly_k300` | 0.000 | **0.000** | 2.889 | 0.030 | 0.028 | **0.894** | **0.092** | **0.101** |
| `blend75` | 0.016 | **0.694** | 2.804 | 0.103 | 0.042 | 0.980 | 0.006 | 0.001 |

**AMBER's earned role reproduces for the sixth time**: the input carries 1.113 close contacts under
2.6 Å and `bond_strain` 0.131; k = 30 takes both to **0.000** and **0.011**.

## 4.4 THE PRE-REGISTERED TEST, AND WHY IT FIRES

The bar (`PREREG_C.md` §4.4): *validity indistinguishable from or better than k30's, **and** ΔCα
smaller than k30's by more than **that comparison's own MDE** (§1.2; the brief's blanket 0.084 Å is
not a valid threshold), fold CI excluding zero.*

**SIX arms beat k30 on Cα past their own per-comparison MDE, and every one of them fails on
validity with a CI excluding zero.** (Under the brief's blanket 0.084 Å constant only three would
have qualified; the correction in §1.2 adds `k100`, `k300` and `caonly_k100` and makes the
conclusion stronger, not weaker.)

| arm | ΔCα vs k30 | SE | own MDE | what disqualifies it |
|---|---|---|---|---|
| `blend75` | −0.1193 | 0.0118 | 0.033 | **clash 2.6 Å +11.63 [+9.67, +14.51]**, bond_strain +0.229 |
| `caonly_k300` | −0.1100 | 0.0109 | 0.030 | **cis_frac +0.4242 [+0.3799, +0.4732]**, rama_fav −0.187, rama_out +0.202 |
| `blend50` | −0.0922 | 0.0082 | 0.023 | **clash 2.6 Å +7.49 [+6.13, +9.34]**, bond_strain +0.156 |
| `caonly_k100` | −0.0776 | 0.0089 | 0.025 | **cis_frac +0.3152 [+0.2695, +0.3476]**, rama_out +0.138 |
| `k300` | −0.0737 | 0.0087 | 0.024 | **clash 2.6 Å +0.971 [+0.784, +1.090]**, bond_strain +0.052 |
| `k100` | −0.0386 | 0.0057 | 0.016 | **clash 2.6 Å +0.561 [+0.471, +0.659]**, bond_strain +0.024 |

**Not one arm buys Cα accuracy for free.**

`blend50/75` are not minimisations at all — they interpolate back toward the unrepaired input, so
they recover the input's Cα by recovering the input's clashes. `caonly_k300` is §4.5.

**And in the other direction, the three arms that are genuinely indistinguishable from k30 on
accuracy are indistinguishable at their own resolution too**: `bbo_k30` +0.0106 (MDE 0.019),
`stage_300_30` +0.0078 (MDE 0.012), `caonly_k30` −0.0054 (MDE 0.017) — all below their own MDE
**and** below the 0.0142 Å rotated-frame floor (§4.6), so they are **NOT MEASURED**, not "matched".
`k10` is the one I previously mis-filed here: at +0.0238 with MDE 0.015 it is **above** its own
resolution and is a genuine, measured degradation.

> **F-C2 fires. The frontier is a single monotone trade between displacement and repair; the
> incumbent k = 30 sits on it, not below it.** This reproduces Sprint 18's *"the repair tax is
> irreducible"* by a completely different route — Sprint 18 tested four native-free spacing
> corrections, this tests **sixteen constrained minimisation protocols spanning three orders of
> magnitude of restraint strength and three restrained atom sets** — and lands in the same place.
> **AMBER's standing role is unchanged: stereochemical repair only, at a displacement cost.**

## 4.5 THE RESULT I WOULD HAVE PUBLISHED IF THE BRIEF ALLOWED A SCALAR VALIDITY SCORE

`caonly_k300` — Cα pinned hard, N/C/O and side chains free — is:

* **0.110 Å [0.092, 0.136] more accurate than the incumbent**, 95W/28L, n = 123, comfortably past
  its own MDE (0.030 Å, SE 0.0109) by 3.6× and the rotated-frame numerical floor by 8×; **and**
* **statistically indistinguishable from the incumbent on clash count**: −0.008 [−0.050, +0.032]
  contacts under 2.6 Å, both effectively zero; **and**
* better than the *input* on every steric measure, with `min_heavy` 2.889 vs 2.752.

**A scalar validity score built on clashes — the most natural single number, and the one this
programme's earlier repair benchmarks leaned on — would have promoted it as the sprint's headline.**

What kills it is elsewhere in the vector: **`cis_frac` +0.4242 [+0.3799, +0.4732]** — 42% more cis
peptide bonds — with `rama_favoured` −0.187 and `rama_outlier` +0.202. The mechanism is immediate
once seen: restraining **only** Cα leaves N and C free, so the peptide plane can flip. And it is
**dose-dependent in exactly the way that confirms the mechanism**: as the Cα pin hardens, cis_frac
runs **0.010 → 0.051 → 0.101** at k = 30 → 100 → 300, while the clash count stays at zero
throughout. The structure has no clashes because it is compact and well spaced. It is also not a
protein.

> **This is the sprint's cleanest vindication of §6 of the brief — report validity as a vector,
> never as a scalar.** One arm, one omitted axis, and the conclusion inverts.

`bbo_k30` is the mirror image and worth recording: restraining N, CA, C **and O** gives the best
backbone dihedrals in the table (`rama_favoured` +0.0467 [+0.0073, +0.0812] and `cis_frac` −0.0522
[−0.0601, −0.0428] against k30) at no measurable accuracy cost — but **worse sterics**
(clash 2.6 Å +0.4390 [+0.3554, +0.4959]). It trades one validity axis for another. Neither arm
dominates; that is the whole shape of this frontier.

## 4.6 THE ROTATED-FRAME NULL, reported with its MAXIMUM as §7 of the brief requires

    n = 15    MAX |delta Ca-RMSD| = 0.014192 A    mean 0.001726 A    max |delta E| 4.769 kcal/mol

Exactly zero by construction (ff14SB, the positional restraint, the projection and `ca_rmsd` are all
rigid-invariant), so this is the pipeline's numerical floor. **Every effect quoted above 0.024 Å is
above the floor; the three "ties" at 0.005–0.011 Å are below it and are reported as NOT MEASURED
rather than as matches.** The 4.77 kcal/mol energy excursion is the known
minimiser-under-translation sensitivity (`verify/amber_audit.json`) and is why the convergence gate
sits at 1000 kcal/mol rather than at a tight threshold.

## 4.7 No iteration bound is configured in this run

`steps = 0` — the deployed unbounded protocol — so there is nothing to fire and nothing to certify.
The bound was **withdrawn**; see §4.8.

## 4.8 A CLAIM OF MY OWN THAT I PUT INTO THE PROGRAMME'S LEDGER AND THEN REFUTED

> ### RETRACTED: *"`core.amber`'s minimisation does not terminate on some inputs."* My own claim, my own error, refuted by my own retest.

**What I claimed.** That `core.amber` calls `LocalEnergyMinimizer` with `steps = 0` (unbounded) and
that on some inputs it never returns — ">40 minutes of a full core at k = 3 on 1A1P", ">35 minutes
on 1D0W", ">7 minutes for a *capped* k = 10 call on 1CEK" — and that it was the **softness of the
restraint** that triggered it. I recorded this as a defect in shared machinery affecting five
sprints, added a uniform iteration bound, and certified the bound with a new gate.

**What is true.** Re-measured on a quiet box, uncapped, on the identical targets and arms:

    1A1P  k=1   wall 14.7 s  cpu 14.2 s  converged
    1A1P  k=3   wall 15.2 s  cpu 14.4 s  converged
    1A1P  k=10  wall 12.8 s  cpu 11.8 s  converged
    1CEK  k=1/3/10        ~6 s each      converged
    1D0W  k=1/3/10   8.7 / 10.0 / 7.1 s  converged

**There is no pathology and no non-termination.** What I measured was **CPU starvation** — six other
lanes' processes on a shared eight-core box — and I attributed a scheduling artefact to the
minimiser. The completed run confirms it: over ~2000 minimisations the **mean per-call wall is
6.5 s and the p95 is 13.6 s**, including every arm I had called pathological.

**How it was caught.** Not by review — by the contradiction between my own claim and my own
completed table, and then by running the direct test rather than arguing the discrepancy away.

**What follows, and it is not free.**

1. **Three pre-registration deviations are WITHDRAWN**, because every one rested on the refuted
   claim: the ladder is restored to the full pre-registered `(1, 3, 10, 30, 100, 300, 1000)`,
   `caonly` to its full declared `(30, 100, 300)`, and the iteration bound is removed. The run above
   is the pre-registration as written.
2. **The capped run was discarded, not patched.** It is preserved as
   `_SUPERSEDED_agentC_pareto_capped2000_n126.json` — complete, valid data — and is quoted nowhere.
   The reason it could not simply be reported is §4.9.
3. **The Z6 rule survives** and is now on firmer ground than the claim that produced it.

## 4.9 THE SECOND ERROR, WHICH THE FIRST ONE HID — a bound must be checked on the calls it BOUND

Gate GC3 certified the `steps = 2000` bound as bit-inert at **0.00e+00 Å** on five targets. Those
five were targets where the incumbent **converges before the bound**, so the certificate was about
calls the bound never touched. The single call in the whole run that the wall-clock proxy flagged as
cap-limited (`k30` on `8T61`) was never checked by the gate. Checked afterwards, directly:

    8T61  k30  capped(2000) vs uncapped:  max delta Ca = 1.595e-01 A,  delta E = 0.418 kcal/mol

**The bound bound, and it moved the structure by 0.16 Å.** One row of a 126-row table was a
different operator from the incumbent, and the gate that was supposed to catch exactly that had
certified the opposite.

Arguing that 0.16 Å on 1 of 123 targets shifts a mean by ~0.001 Å would have been true, and would
have been the wrong call: **the confound was removable, so it was removed** by rerunning unbounded.

> **The rule, which generalises past this instance and past Z6:** an inertness certificate must be
> evaluated **on the calls where the intervention actually acted**. A certificate computed only on
> the cases the intervention did not touch is not weak evidence — it is *no* evidence, and it will
> read as a pass.

---

# 5. DEVIATIONS FROM THE PRE-REGISTRATION — every one, with its reason and its timing

| # | declared | done | why, and when the decision was taken |
|---|---|---|---|
| 1 | primary condition (1): `mean Δ(readout²) > 0` with a CI excluding zero | reported as **failing as written** (+0.4671 [−0.0308, +1.1057]) beside the metric-scale damage (+0.0631 [+0.0097, +0.1236]) | a mis-specification of my own: I put the primary on the squared scale so the channel shares would be exact, and that scale is too heavy-tailed to clear at n = 126. **Kept, not edited.** Both scales printed. |
| 2 | §4.3 wider restrained set `{N, CA, C, O, CB}` | run as **`bbo` = {N, CA, C, O}** | **CB does not exist for glycine** and `core.amber._index_topology` builds `_restraint_idx` by direct dict lookup, so the declared set raises `KeyError` on most targets. Renamed so no table can read as the declared arm. Decided at implementation time, before any number. |
| 3 | §4.3 restraint ladder `(1, 3, 10, 30, 100, 300, 1000)` | **WITHDRAWN — the full declared ladder was run** | I first dropped k = 1 and k = 3, citing calls that "had not returned after 40 minutes". **That timing fact was CPU starvation, not the minimiser** (§4.8): re-measured on a quiet box the same calls take 14.7 s and 15.2 s. The deviation is void and both rungs were run. They turn out to be **worse** than the incumbent (k = 1 is +0.090 Å with more Ramachandran outliers), so the deviation cost the sprint nothing scientifically — but it was still made on a false premise, and it is recorded as one. |
| 4 | §4.3 `caonly` at k = 30/100/300 | **WITHDRAWN — all three were run** | rested on #3. |
| 5 | rotated-frame null on all targets (n unspecified) | run on the **first 15 targets**, reported with its **maximum** | it is exactly zero by construction, so it is a numerical-floor measurement, not an effect estimate. |
| 5b | `steps = 0` (the incumbent's unbounded minimisation) on every arm | **WITHDRAWN — the deployed `steps = 0` is used on every arm, no bound anywhere** | I introduced a uniform iteration bound to fix a non-termination defect **that does not exist** (§4.8), and the bound then perturbed one measured call by **0.16 Å** (§4.9). Both the bound and the deviation are removed; the capped run is preserved as `_SUPERSEDED_*` and quoted nowhere. |
| 6 | §3 arms only | **two declared extensions added** (`agentC_kv2`, `agentC_kv3`) | each is labelled in its own module docstring with the result it was written after, and each could only make the mechanism claim harder to sustain. §3.4's primary was not re-run, re-cut or re-scaled. |

**Deviations 3, 4 and 5b were all WITHDRAWN once the timing claim underneath them was refuted, and
the final Q2 run is the pre-registration exactly as written.** Nothing was re-normalised after an
unfavourable answer, no ladder was extended past its pre-registration, no sign was flipped post hoc,
no easier metric was substituted, and the sealed 60-target benchmark was not read, probed, derived
from or tuned against.

---

# 6. LEDGER

| # | hypothesis | mechanism | primary metric | controls | n | result | status |
|---|---|---|---|---|---|---|---|
| C1 | a score gate damages the average by **collapsing ensemble ambiguity** | Krogh–Vedelsby `readout² = E_mem² − D²` | Δ(readout²), −Δ(D²), share `S_D` | matched-random ×30; zero-information helix; **score-free `rand_lowD` at matched D** | 126 | −Δ(D²) +0.692 [+0.581, +0.843] and `S_D` = 1.48, **but** a score-free gate at the same D costs +0.003 against Legacy's +0.063 | **REFUTED as a cause** (the correlation is ESTABLISHED) |
| C2 | a **diversity-preserving** Legacy gate recovers the deficit | keep each Cα-RMSD cluster's best-Legacy member | Cα-RMSD of the emitted average | plain `legacy`, matched-random, no gate | 126 | −0.0664 [−0.1392, −0.0087] vs `legacy`; −0.0033 [−0.0198, +0.0112] vs random | **SUPPORTED vs Legacy; NOT MEASURED vs random** |
| C3 | the damage is the **displacement of the emitted mean along the pool's own error direction** | `readout² = ‖b‖² + 2‖b‖·ALIGN + ‖Δ‖²` | ALIGN excess over a 200-draw matched-random null | four score-free gates spanning the whole diversity range | 126 | score gates +0.0209 [+0.0004, +0.0402]; score-free −0.0048 [−0.0124, +0.0020]; difference −0.0256 [−0.0455, −0.0008]; per-target r = +0.95 | **SUPPORTED** |
| C4 | `‖Δ‖` is a **native-free design-time proxy** for gate damage | the magnitude term of C3 | Spearman over gate designs | — | 12 designs × 126 | ρ = **+0.930** across designs; r = +0.169 per (arm, target) | **SUPPORTED at design level, weak per target** |
| C5 | Legacy `steric` can reject **impossible** candidates | hard steric rejection of the worst `r` | Cα-RMSD of the emitted average, r = 10 | matched-random ×30 **and** a diversity-preserving rejection | 126 | +0.0175 [+0.0048, +0.0348] vs random, +0.0170 [+0.0091, +0.0286] vs `divkeep`, +0.0204 [+0.0081, +0.0376] vs no rejection | **REFUTED, with a sign. CLOSED.** |
| C6 | a differently-constrained AMBER buys the same repair for less Cα displacement | restraint strength, restrained atom set, staging, pull-back, projection | Cα-RMSD **and** the validity vector, as two axes | its own gated input; the incumbent k30; a rotated-frame null with its MAXIMUM | 126 | every arm matching k30's validity is within ±0.024 Å of it; the three beating it past the MDE carry +11.63/+7.49 clashes or +0.424 cis | **REFUTED — F-C2 fires** (§4.4) |
| C6b | a **scalar** validity score suffices | clash count as the single validity number | Cα and the full validity vector | the vector itself | 126 | `caonly_k300` −0.110 Å with clash indistinguishable from k30, disqualified only by `cis_frac` +0.4242, dose-dependent in the Cα pin | **REFUTED — vindicates the vector rule** (§4.5) |
| C7 | `core.amber`'s minimisation does not terminate on some inputs | soft restraint + unbounded `steps = 0` | wall AND cpu clock of one call, on a quiet box | the same calls uncapped when the box is idle | 3 targets × 3 k | **14.7 / 15.2 / 12.8 s, all converged** — the original >40 min was CPU starvation from six other lanes | **RETRACTED — my own claim, my own error** (§4.8) |
| C7b | the brief's instrument-wide 0.084 Å MDE is a valid threshold | `MDE = 2.8016·SE`, SE from the paired difference | per-comparison MDE for every decisive contrast | — | 126 | per-comparison MDE spans **0.010–0.113 Å**; the constant hid 4 measured Q2 arms and understated a Q1 null's resolution; **no verdict reverses** | **REFUTED — MDE is per-comparison** (§1.2) |
| C8 | an inertness certificate is evidence wherever it is computed | GC3 certified the bound on converging calls only | bit-identity on the call the bound actually bound | the flagged call, checked directly | 1 | **0.1595 Å / 0.418 kcal-mol** — the bound bound, and the gate had certified the opposite | **REFUTED — a certificate must be evaluated where the intervention acted** (§4.9) |

---

# 7. WHAT IS OPEN, HONESTLY

1. **Does the ALIGN mechanism survive the terminal AMBER repair?** Still the top open item, and now
   the only one from Q2's side too. Everything in §2 is measured on
   the coordinate average. Sprint 18 showed the filter ordering survives an AMBER pass
   (`legacy`→AMBER stays +0.078 worse than `none`→AMBER), so it probably does — but
   `legacy_clust` → projection → AMBER was **not** run here, so the conversion is
   **NOT MEASURED**. It is ~40 minutes of AMBER on three arms and it is the obvious next check.
   This programme has been burned before by a gate whose gain did not convert (`leg_torsion`).
2. **Why is a physics score correlated with the pool's error direction at all?** §2.5 measures that
   it is, not why. `helix` being the worst arm in the table says the answer is probably not about
   physics: any rule that pulls the survivor set toward a single compact reference structure appears
   to move the emitted mean the way the pool is already wrong. Whether that is a property of the
   retrieval pool or of the objective is untested.
3. ~~The soft end of the restraint ladder is untested.~~ **CLOSED — it was tested and it is worse.**
   `k = 1` is **+0.090 Å [+0.057, +0.129]** worse than the incumbent with **more** Ramachandran
   outliers (+0.0340 [+0.0170, +0.0517]). Softening below k = 30 costs accuracy and buys no
   validity. There is no better Pareto point at the soft end.
4. **`leg_contact` is a different animal from the other four gates** and this sprint only noted it:
   it is the one score whose damage is member quality rather than displacement, consistent with
   Sprint 17/18's anti-ranking measurement and with nothing else measured here.

---

# 8. REPRODUCTION

    python -m s19.agentC_lib                # gates GC0, GC1                        ~1 min
    python -m s19.agentC_pareto --gate      # gate GC2                              ~1 min
    python -m s19.agentC_kv                 # Q1, the pre-registered primary        ~2 min
    python -m s19.agentC_kv2                # Q1, declared extension 1              ~9 min
    python -m s19.agentC_kv3                # Q1, declared extension 2              ~7 min
    python -m s19.agentC_reject             # Q3                                    ~12 min
    python -m s19.agentC_pareto --gate      # gates GC2, GC3                        ~5 min
    python -m s19.agentC_pareto             # Q2, n = 126, 16 arms, steps = 0       ~4.3 h
    python -m s19.agentC_report  kv         # every Q1 table
    python -m s19.agentC_report  reject     # every Q3 table
    python -m s19.agentC_report  pareto     # every Q2 table
    python -m s19.agentC_kv2 --report ; python -m s19.agentC_kv3 --report

Every artefact carries `complete`, a row count and a config hash; `s18.phys_lib.read_complete`
refuses to read a partial as a result. Smokes are written to `_SMOKE_*` and partials to
`_PARTIAL_*`, and neither is quoted anywhere above. `OMP_NUM_THREADS=1` throughout.

## Instrument integrity

| check | result |
|---|---|
| sealed benchmark | **not read, probed, derived from or tuned against** |
| metric | full-chain Cα-RMSD, frozen `s12.instrument.ca_rmsd`, all residues, proper rotations, model 1; no alternative reported as primary |
| native information | evaluation and labelled ORACLE diagnostics only. `D`, `‖Δ‖` and **every gate and every rejection rule** are native-free; `readout`, `E_mem`, `FRAME`, `ALIGN` are ORACLE |
| seeding | `s15.seed.stable_rng` only; no bare `hash()`, no unseeded `np.random` in any Sprint-19 Agent-C module |
| Legacy weights | `DEFAULT_WEIGHTS`, never fitted |
| AMBER discipline | convergence gate declared in `PREREG_C.md` §2 before use and reported with its exclusion count; every gated arm compared to **its own gated input**; rotated-frame null reported with its **maximum** |
| cross-sprint reproduction | Legacy components and ORACLE labels identical to `s18/results/down.json` at **0.00e+00**; the Sprint-18 gate table reproduces to the third decimal on all five scores with a 10× tighter null; the operator-law misses reproduce to the third decimal |
| statistics | target-level paired bootstrap, **fold-clustered CI quoted as primary**, median and W/L beside every mean, and — after the correction in §1.2 — **SE and the per-comparison MDE beside every decisive mean**. The brief's instrument-wide 0.084 Å constant is not used as a threshold anywhere |
| completion flags | every quoted artefact `complete` |
