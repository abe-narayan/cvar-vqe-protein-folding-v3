# WORKSTREAM D — THEORY, LITERATURE AND ADVERSARIAL AUDIT
## Sprint 21 findings

Pre-registered in `s21/PREREG_D.md`, written before the corresponding results were inspected and
**not edited after**. Labels per `s21/BRIEF.md` §9. **TARGET is the unit** throughout. MDE at 80%
power = **0.084 Å**. Every structural number states its basis — *point cloud* / *built chain* /
*repaired emission*.

**Benchmark seal, proved by hash and not by read, before anything else was done:**

    results/benchmark_manifest.json
      sha256  a40581ad01cfd2b77aa755a149fe40f3661dc5d35d4c4be0133f6bf23248422d
      size    8002 bytes      mtime 2026-09-01 18:52:50

Byte-identical to Sprint 20's certificate. The file was not opened, probed or derived from.

**Where each pre-registration lives.** D1–D4 are in `s21/PREREG_D.md`, written before any of their
results were inspected and **not edited after**. D7, D13, D14, R1 and the `latentfull` audit were
registered in their own module docstrings — hypothesis, prediction, primary endpoint, falsifier,
null, matched control, budget — before those modules were run, which is the same discipline applied
to work that arose during the sprint. **Two of my six registered hypotheses died on their own
falsifiers** (R1 and D7) and are reported as such.

Artefacts, each with a completion flag that requires the **full** configuration:

    d_cvarop.json          .COMPLETE   the readout operator, n = 126, 3 H x 3 alpha x 4 readouts
    d_sim.json             .COMPLETE   exact enumeration of the deployed ansatz, 64 cells
    d_null.json            .COMPLETE   the two matched-displacement nulls, n = 30
    d_divsel.json          .COMPLETE   R1, set-diversity selection, n = 126 -- my own hypothesis REFUTED
    d_mde.json             .COMPLETE   the MDE derived per comparison, 26 comparisons
    d_partial.json         .COMPLETE   Legacy's rank skill partialled on the distogram, n = 42
    d_distobj.json         .COMPLETE   the two 'distance' functionals, n = 126 -- my own hypothesis REFUTED
    d_lrank.json           .COMPLETE   native-free readouts inside the objective's top-M, n = 75
    d_enc_gate.json        .COMPLETE   48/48 bit-for-bit vs Sprint 20's own stored rows
    d_enc_LEG.json         .COMPLETE   encoding x budget grid, LEG, n = 10, 4 seeds
    d_enc_AMBc_spsa_B512   .COMPLETE   the decisive encoding cell, n = 10, 4 seeds, 80 cells
    d_enc_LEG-DIST_n60_... .COMPLETE   the POWERED encoding replication, n = 60, 4 seeds
    _PARTIAL_d_enc_AMBc_fullgrid_3targets.json   the abandoned full grid, preserved not deleted

All under `s21/results/`. Code: `s21/d_cvarop.py`, `d_enc.py`, `d_sim.py`, `d_null.py`,
`d_divsel.py`, `d_mde.py`, `d_partial.py`, `d_distobj.py`, `d_lrank.py`.

---

# 0. EXECUTIVE SUMMARY

> ## The sprint's opening endpoint is maximised by an energy containing no information; the pool-restricted tail-average selector *is* the incumbent pipeline, bit-exact; and no quantum-resource claim is reachable on a 9–16 qubit register at any bond dimension or depth.

| # | statement | label |
|---|---|---|
| **D1** | **`tail_avg(α=0.15)` under the distogram is the shipped pipeline, bit-exact.** 3.0483 Å over 126 targets, `max |diff| = 0.000000` against `bench_results/cache/1fc9f2dcf489e2fb`'s `rmsd_avg` on **126/126** targets, corr 1.000000. And `argmin\|disto` = 3.454 reproduces `s12/instrument`'s asserted `shipped argmin 3.4540`. The "CVaR tail vs argmin" question was settled in production five sprints ago: the pipeline already takes the tail-average. | **EXACT** (identity), verified n = 126 |
| **D2** | **The pre-registered primary endpoint `tail_avg(0.15) − argmin` is maximised by a ZERO-INFORMATION energy.** `rand` (a per-target random permutation, 16 independent draws) gives **−1.035 [−1.126, −0.946]**, larger than `legacy` (−0.735) and `disto` (−0.406). 181% of the physics mean; 1.41× the weakest physics arm. My falsifier does not fire. The endpoint measures the readout, not the Hamiltonian. | **ESTABLISHED**, n = 126 |
| **D3** | **A Hamiltonian's sign flips with the readout.** Against matched-count random controls, Legacy **beats** them on `tail_member` (−0.41 [−0.56, −0.26]) and **loses** to them on `tail_medoid` (+0.30) and `tail_avg` (+0.33), all CIs excluding zero. "Which Hamiltonian is best" is not well defined without fixing the readout. | **ESTABLISHED**, n = 126 |
| **D4** | **Pool-restricted CVaR with an order-based readout cannot beat argmin.** `tail_min ≡ argmin` in 9/9 (H, α) cells, max difference 0. Q6 puts the CVaR-optimal law on the α-tail; the α-tail contains the pool argmin; the deployed readouts (`core/quantum.py:1600–1626`) are all order-based. Scope: pool restriction, same H for training and readout, order-based readout. | **EXACT** (with scope) |
| **D15** | **At Sprint 20's 8192-evaluation budget, the budget equals or exceeds the ENTIRE latent space on 75/126 targets** (one qubit per residue, two basins each, so \|latent\| = 2^n; median 8192). On a majority of targets the VQE is resampling a space it could have enumerated within budget, so Q3 is **not** a search failure there — it is an ordering failure, now with an exhaustiveness argument rather than only six instruments. | **ESTABLISHED** (arithmetic + source) |
| **D5** | **No quantum-resource claim is reachable here, for any ansatz.** The deployed circuit is `MPSAnsatz(n, layers=2, ...)` with **n = residues**, so the register is **9–16 qubits**. I enumerated the circuit's complete law exactly at every register size and at layers 1…20, including χ saturated at χ_max: **n = 16, layers = 20, 65,536 states, 64 ms, Σp = 1.0000000000.** χ = 4 is true and is not the binding constraint. | **ESTABLISHED by measurement** |
| **D6** | **Two of Sprint 20's "twelve encoding cells" are one measurement.** `AMB\|nelder` ≡ `AMBc\|nelder` bit-identical on 10 targets × 2 seeds in both encodings — `AMBc` is a monotone transform and Nelder–Mead is comparison-only. Sign test p = 0.0063 → 0.0117 on 11 distinct cells, and 2 (not 3) distinct cells have CIs excluding zero. | **EXACT** (the identity) / labelling correction |
| **D7** | **The encoding lever is confounded with step count in every cell,** and for SPSA the confound is a *diagnostic*: `track_quality`'s FD probes cost `2·d` and `d` doubles in the embedding, eating **288/512** budget units embedded against 144/512 in θ. With the diagnostic off, SPSA's iteration counts are **exactly equal** (256 vs 256). | **ESTABLISHED** (from source + my own run) |
| **D8** | **The encoding lever does not survive the step-count control.** `AMBc/spsa`, the largest Sprint-20 effect and one of only two distinct significant cells, goes from **−0.649 [−1.205, −0.165] "sig"** to **−0.295 [−0.636, +0.022]**, 6W/4L, at *exactly equal* step counts (256/256) with 4 seeds. On the **powered replication (n = 60, 3840 cells, COMPLETE)** the pooled effect is **−0.0061 [−0.0448, +0.0329], W/L 30/30, sign p = 1.0000**, own MDE ≈ 0.055 — a **powered** flat null. Where steps differ the "significant" effect **flips sign between objectives** (LEG `adam` −0.130, DIST `adam` **+0.091**, both CIs excluding zero) and both collapse under the step control. | **NOT SUPPORTED** |
| **D19** | **The encoding effect is quantitatively PREDICTED by the objective's own step-count response,** with no free parameter: predicted −0.0912 vs observed −0.1301 on LEG, and predicted **+0.0865 vs observed +0.0905** on DIST — residual **+0.0040 [−0.023, +0.027]**. Both residuals span zero. **And this refutes my own stated mechanism**: the step response is significantly signed in **both** directions (LEG +0.0912, DIST −0.0865, both CIs excluding zero), so *"optimising harder hurts"* is a property of the **Legacy** objective, not of the instrument — exactly as `s21/BRIEF.md` §5 already states and I failed to. | **ESTABLISHED** / my own mechanism **REFUTED** |
| **D16** | **`_EmbField` has a gauge bug, and it is the durable output of D1.** `‖u‖` grows from 1.000 to **2.50 (max 4.03)** along a direction whose true gradient is **exactly zero**, so the embedded arm anneals its own step size by ~2× and moves **10% less far** in torsion space at identical step count. An embedded arm and a θ arm are running different **optimisers**, not different **coordinates**. Fix: project the update onto the tangent, or renormalise `u` each iteration. | **ESTABLISHED** |
| **D9** | **The gauge radius is a real, large, unintended annealing schedule.** `‖u‖` starts at 1.000 and reaches **1.94 (max 3.06)** for SPSA@512 and **2.67 (max 3.79)** at B = 256, in a direction whose true gradient is exactly 0 because `arctan2` is scale-invariant. The effective angular step decays as 1/‖u‖. | **ESTABLISHED** |
| **D10** | **"72% of AMBER's damage is move size" is a ratio against a null that itself carries 0.554 Å of skill.** The two matched-magnitude nulls differ by **−0.5540 [−0.8310, −0.2879] fold-clustered, 24W/6L**. Against the isotropic null the same arithmetic reads **161%** and AMBER's direction is **better** than noise by −0.378. The number is right; the label "zero information" on `toward_member` is wrong. | **ESTABLISHED** — labelling correction |
| **D11** | `s21/tailprice.py` as first committed **never ran** (two wrong signatures; one of them reaches for the discrete lattice) and its artefact was stamped `complete: true` with **3/3 rows skipped**. Corrected by the coordinator within the hour. | recorded as process |
| **D14** | **Legacy carries NO in-band rank information beyond the distogram.** Its marginal Spearman against true RMSD is **+0.3281 [+0.3011, +0.3552]** fold-clustered; partialled on the distogram — with which it correlates **+0.4784** — it is **−0.0076 [−0.0616, +0.0681]**, median −0.033, 25/42 targets negative. It retains **−2%**. Control: AMBER's partial stays ~0, so the partialling is not manufacturing structure. | **ESTABLISHED** (n = 42, Workstream A's live artefact) |
| **D13** | **"MDE = 0.084 Å" is not a property of the instrument.** MDE = 2.8016·sd(paired differences)/√n, so it belongs to the *comparison*. Across 26 real comparisons the per-comparison MDE spans **0.11× to 9.54×** the quoted constant — a factor of **84** — and only 12% fall within ±1.5× of it. It errs both ways: the deployed AMBER repair has MDE **0.0096 Å** and a **6.09 SE** effect the briefs call "a quarter of the MDE"; Sprint 20's `AMBc/spsa` encoding cell has MDE **0.8014 Å**, larger than the −0.649 it declared significant (power 0.62, Type-M 1.27×). | **ESTABLISHED** — methodological |
| **D18** | **The exhaustive-latent control (§3.1a's O7) is built, run and COMPLETE at n = 126 — it closes the search half by construction, on the whole instrument.** Given exhaustive access to **every one** of `2^n` latent configurations, the objective's **exact** argmin still misses the latent's own ORACLE best by **+1.8149 [+1.6129, +2.0594] with 0 of 122 targets reversed** (+1.7202, 0/119, under the Bayes score). The latent contains structures **1.18 Å better than the incumbent** on 104/126, and the objective **does** order it (beats the space's own mean by −1.3674 on 109/126). The matched primary is **+0.0851 [−0.0176, +0.1691]**, 53W/73L — magnitude **NOT MEASURED**, directional falsifier **DID NOT FIRE**. | **ESTABLISHED**, n = 126 |
| **D25** | **The objective's argmin over 512 RANDOM latent draws is statistically indistinguishable from its argmin over the ENTIRE 8192-configuration latent** — `rand{M}_argmin` falls 4.665 → 4.040 → 3.684 → **3.494** against the exhaustive 3.449, ending at **−0.045 [−0.177, +0.137], NOT MEASURED.** A third independent route to search saturation, and the most operationally direct: the search is *done* at 512 draws, measured against its own exact optimum. **The difference between 512 and 8192 evaluations is under the MDE.** | **ESTABLISHED**, n = 45 interim → final in the artefact |
| **D26** | **T2 is REFUTED, and my error-coherence mechanism with it.** The latent's top-M averages against its matched random control by **−0.478** (M = 64) and **−0.327** (M = 512) against the pool's **−0.377**, clearing its own MDE at every M — comparable or larger, not smaller. **The T2 generalisation I offered for the next sprint is WITHDRAWN** (§4.35) rather than rescued with an unregistered selectivity caveat. T1 and T3 SUPPORTED. | **REFUTED** (my hypothesis) |
| **D22** | **The source swap loses, and the fork verdicts I gave were computed rather than requested.** `lat_avg75 − pool_avg75` (basis-matched) = **+0.3898 [+0.3225, +0.4793]**, 47W/79L at n = 126: the latent through the *identical shipped operator* is 0.41 Å worse than the pool, and its ORACLE ceiling is 0.56 Å worse too. The composite footnote decomposes as **96.2% source, 3.8% basis** (basis price +0.0155 [+0.0093, +0.0232]). The `n = 9` hazard does not move the primary — which instead **grows monotonically** with the length floor (+0.427 at n ≥ 10, +0.475 at n ≥ 12), the min-of-M lesson in a third place. | **ESTABLISHED**, n = 126 |
| **D23** | **My consistency prediction is REFUTED and the reconciliation is the finding.** `lat_avg75 − lat_avg75_rand = −0.4290 [−0.5077, −0.3317]`, 91W/35L, n = 126. D9 compares the set's **BEST** (no skill past M ≈ 8); `latentsel` compares the set's **MEAN** through averaging (0.43 Å of skill). **The objective's ordering is useless for finding the best member and worth 0.43 Å for shaping the mean — which matters is decided by the READOUT.** My §1.2d, third arrival. | **ESTABLISHED** / my prediction **REFUTED** |
| **D24** | **In L19's bimodal decomposition the FAILURE cell's sign is forced by the split.** Conditioning on "the top-512 misses the ORACLE-best" selects targets where that arm failed and imposes no matching condition on the random arm. The **magnitude** survives (top-512 sits 1.23 Å above the latent ORACLE where a random window sits 0.38 Å above) but the **sign is guaranteed**. The clean version is unconditional excess-over-ORACLE for both windows. | **ESTABLISHED** — methodological |
| **D21** | **The top-M ceiling ladder was min-of-M, and the control I registered before its data existed overturned the headline.** `GAP = rand_M − top_M` is **−1.184 [−1.448, −0.918] at M = 1** and **+0.107 [−0.051, +0.300] at M = 512**: the objective's ordering skill is concentrated at M = 1 and exhausted by M ≈ 8. 512 **random** latent draws reach **1.879** against the objective-ordered 512's 1.986. **The objective is a good FILTER and a useless RANKER.** My own prediction that the gap would be substantial at every M is **REFUTED**. | **ESTABLISHED** (coordinator's measurement, my registered control) / my prediction **REFUTED** |
| **D20** | **Budget-versus-latent-size was not the operative variable.** On the 51 targets whose latent *exceeds* the deployed 8192-evaluation budget — the only ones where the sampler was genuinely searching rather than resampling — the matched primary is **+0.1738 [−0.0079, +0.3715]**, if anything *worse* than on the enumerable-within-budget half. My P2 and P3 confirmed; **my P1's magnitude is NOT MEASURED and its mechanism REFUTED** (`corr(n, latent_oracle) = +0.297`: the ORACLE ceiling gets *worse* with `n` despite a 128× larger set). | **ESTABLISHED** / my mechanism **REFUTED** |
| **D17** | **The brief's §2 opening row is a 20-target subset, and that subset is harder than the instrument.** At n = 126 on the window basis the shipped selector's argmin is **3.4540** (= `s12/instrument`'s own asserted constant) and the pool mean **4.4533**; the brief quotes 3.676 and 4.739. On Sprint 20's own `subset(20)` I reproduce **pool mean 4.7318 and pool best 1.8882** against its reported 4.739 / 1.901 — agreement to 0.01 Å. **The gap is the subset, not the functional and not the basis.** And my own hypothesis that the two distance functionals differ is **REFUTED**: median `ρ(E_bayes, E_sq) = 0.973` and their argmins agree to +0.050 Å against that comparison's own MDE of 0.124. | **REFUTED** (my hypothesis) / provenance correction |
| **D12** | **My own Priority-4 proposal is REFUTED at n = 126, and its sign control inverts.** Selecting the averaging set for maximal geometric diversity is dead on all three score bands; *minimising* diversity **beats** its matched control by −0.214 [−0.397, −0.040]. The mechanism behind D3 is therefore **error coherence, not set spread** — corrected in §1.2d rather than left standing. | **REFUTED** (my own hypothesis) |

---

# 1. PRIORITY 1 — ATTACKS ON THIS SPRINT'S LIVE CLAIMS

## 1.1 The sprint's opening move did not run, and said it had

`s21/results/tailprice.json` at 09:06 held three rows, **all skipped**, each with
`ValueError("could not convert string to float: 'INWKGIAAMAKKLL'")`, and
`"complete": true`.

Three separable defects, reported to the coordinator within the hour and all three acted on.

**(a) Two wrong signatures.** `s16.energy_lib.legacy_components_of_windows` is
`(seq, PHI, PSI, chunk)`; the call passed `(pdb, W, seq, fold)`, so the sequence string landed in
the `PSI` slot. **And the second one is not a typo you can patch**:
`s13.qarch_lib.amber_energies` is `(space: Space, S, ...)` where `Space` is the **discrete k-state
lattice** (`qarch_lib.py:44–64`) and `S` is an integer state vector fed to
`core.amber.single_point(seq, rep, [int(x) for x in S[b]])`. That interface cannot score real
retrieved windows, which carry continuous torsions, and reaching for it regresses to the lattice
`s21/BRIEF.md` §6 forbids. The continuous route is `AmberHamiltonian.energy_from_coords` /
`core.amber.refine_coords(steps<0)`. Separately, `legacy_components_of_windows` needs `PHI/PSI`
and the module only ever loaded `W`, so **neither physics arm had an input**.

**(b) A completion flag that fired on the wrong condition.** `complete = len(rows) == len(tg)`
counted **skipped** rows, so 3/3 total failures wrote `complete: true`; and it did not require
n = 126, so a 3-target smoke satisfied it. This is exactly the hazard `s21/BRIEF.md` §6 names, and
it is the *fifth* completion-flag defect this programme has caught. Every artefact of mine this
sprint requires the full configuration — all targets, all Hamiltonians, all α, all readouts, all
seeds, all budgets, and every value finite — and **deletes** the flag when the condition fails.

**(c) The scientific defect, which survives fixing (a) and (b).** The pre-registered PRIMARY
endpoint was `tail(0.15) − argmin`. `argmin` returns **one** structure; `tail{a}` returns the
**coordinate average of a·m** structures. Those differ by the averaging operator, which project
memory `averaging-space-beats-the-objective` prices at ~1.0 Å while "any objective contributes only
0.171 Å". So the endpoint was near-certain to "confirm" its prediction for a reason that has
nothing to do with the Hamiltonian. §1.2 measures exactly that.

And the tail arm is not the operator the pillar has. From source, `core/quantum.py:1600–1626`, a
run emits `vqe_bitstring` (argmin over the final distribution's samples), `vqe_modal_bitstring`
(the mode) and `best_seen_bitstring` (argmin over everything seen). **None is a tail mean.** CVaR
is the *training objective*; the *readout* is order-based.

## 1.2 D2 — THE READOUT OPERATOR. n = 126, COMPLETE

`s21/d_cvarop.py` → `d_cvarop.json` + `.COMPLETE`. Pool-restricted (K = 500), no optimisation, the
native read only to score. Four readouts on the same α-tail, each beside its **matched-count**
control, for three Hamiltonians one of which contains nothing.

### 1.2a The instrument reproduces both of the programme's pinned constants exactly

| quantity | mine | the programme's | agreement |
|---|---|---|---|
| `tail_avg(0.15) \| disto` | **3.0483** | `rmsd_avg` in `bench_results/cache/1fc9f2dcf489e2fb` = **3.0483** | `max\|diff\|` = **0.000000** on **126/126**, corr 1.000000 |
| `argmin \| disto` | **3.454** | `s12/instrument`'s asserted `shipped argmin 3.4540` | matches |

It has to. The shipped pipeline **is** "rank the K = 500 pool by the Bayes-risk distogram score,
keep the top 75 (= 0.15 × 500), coordinate-average them". So:

> **The pool-restricted α-tail-average selector at α = 0.15 under the distogram IS the incumbent.**
> The mandatory matrix's `disto` row at α = 0.15 is not a new arm — it is the production pipeline
> wearing a CVaR label, and "argmin vs tail" was settled in production: the pipeline already takes
> the tail-average, which is exactly why it sits at 3.048 (point cloud) and not at the 3.454
> argmin. — **EXACT**, verified at n = 126.

This is also the strongest available soundness gate on everything else in this section.

### 1.2b The table. Note the basis column

Mean over 126 targets. `argmin`, `tail_member`, `tail_medoid` return a **real retrieved window**
(built basis, valid geometry). `tail_avg` returns a **POINT CLOUD** — not a buildable backbone
(s20 X1) — and comparing it to a single-member row is a basis mismatch.

| readout | basis | matched control | disto | legacy | **rand** |
|---|---|---|---|---|---|
| `argmin` | window | 4.434 | 3.454 | 4.490 | 4.453 |
| `tail_member` 0.15 | window | 4.434 | **3.530** | 4.026 | 4.446 |
| `tail_medoid` 0.15 | window | 3.730 | 3.282 | 4.031 | 3.715 |
| `tail_avg` 0.15 | **POINT CLOUD** | 3.425 | **3.048** | 3.754 | 3.418 |
| `ORACLE_tail_avg` 0.15 | POINT CLOUD | — | **1.963** | — | — |

ORACLE selector ceilings at α = 0.05 / 0.15 / 0.30: **1.644 / 1.963 / 2.275** (point cloud).

### 1.2c THE PRIMARY: a zero-information energy wins the endpoint

`tail_avg(0.15) − argmin`, paired over 126 targets. **CI construction quoted**: i.i.d. over
targets, with a fold-clustered interval beside it.

| H | mean | CI95 iid | CI95 fold-clustered | W/L |
|---|---|---|---|---|
| `disto` | −0.406 | [−0.536, −0.276] | [−0.513, −0.317] | 95/31 |
| `legacy` | −0.735 | [−1.113, −0.369] | [−1.108, −0.402] | 92/34 |
| **`rand`** | **−1.035** | **[−1.126, −0.946]** | **[−1.207, −0.914]** | **123/3** |

`rand` is a per-target random permutation of the pool, every readout averaged over **16
independent permutations** so its `argmin` is not a one-draw quantity. It carries **no
information**, and it produces the **largest** gap of the three: 181% of the physics mean, 1.41×
the weakest physics arm. My pre-registered falsifier — *"D2b is REFUTED if `rand`'s gap is less
than half the physics gaps"* — **does not fire**.

> **The endpoint `tail − argmin` cannot rank Hamiltonians. It measures the readout.** — ESTABLISHED

Decomposed with the tail held fixed:

| α | H | `tail_avg − tail_medoid` (the AVERAGING operator) | `tail_medoid − argmin` (SET vs MINIMUM) |
|---|---|---|---|
| 0.15 | disto | −0.234 [−0.304, −0.161] | −0.172 [−0.315, −0.029] |
| 0.15 | legacy | −0.277 [−0.367, −0.191] | −0.458 [−0.890, −0.031] |
| 0.15 | **rand** | **−0.296 [−0.395, −0.195]** | **−0.739 [−0.902, −0.577]** |

The averaging half is −0.23 to −0.31 for **all three**, `rand` included — it is energy-independent
by construction. The set-vs-min half is −0.74 for `rand` because argmin-of-noise is a random pool
member and a medoid is not. **Neither half is physics.**

### 1.2d THE FINDING WORTH KEEPING: a Hamiltonian's sign flips with the readout

Each readout minus its **matched-count** control (α = 0.15; the same pattern holds at 0.05 and
0.30 with every CI excluding zero).

| readout | control | disto | legacy |
|---|---|---|---|
| `tail_member` | random pool member | **−0.904 [−1.078, −0.729]** 108W/18L | **−0.407 [−0.557, −0.259]** 82W/44L |
| `tail_medoid` | medoid of a random size-k subset | **−0.447 [−0.671, −0.241]** 83W/43L | **+0.302 [+0.204, +0.404]** 35W/91L |
| `tail_avg` | average of a random size-k subset | **−0.377 [−0.547, −0.215]** 85W/41L | **+0.329 [+0.211, +0.448]** 41W/85L |

**Control validation**: `rand`'s three arms against the same controls are +0.013, −0.015, −0.007,
all spanning zero. The controls are matched.

> **Legacy's tail is BETTER than a random subset if you draw one member from it, and WORSE than a
> random subset if you take its medoid or its average.**

**The mechanism, corrected by my own R1 falsification (§4.1) — the first version of this sentence
was wrong and is retracted here rather than left standing.** I originally wrote that Legacy's tail
is a *low-diversity* set and that low diversity destroys what the averaging operator lives on. R1
tested exactly that and killed it: **minimising** set diversity **helps** the averaging readout
(−0.214 [−0.397, −0.040] against its matched control). What matters is not the set's **spread** but
whether its errors are **coherent**. Averaging cancels i.i.d. error and preserves systematic error
— `error-coherence-decides-correctors`, where at identical 0.688 sign accuracy coherent mistakes
emit +0.31 Å and i.i.d. mistakes −0.14 Å. Legacy's tail shares a **coherent compactness bias**
(L2c: Legacy-preferred candidates are 0.45 Å more compact, 124W/2L), so averaging preserves the
bias; a random subset's errors are nearer independent, so averaging cancels them. The measurement
in the table stands exactly as reported; only the mechanism moved.

**Consequence for §4 of the brief.** The mandatory matrix must state its readout on every row. The
same matrix run at `tail_avg` reports Legacy as harmful and run at `tail_member` reports it as
helpful. **Both are true.** "Which Hamiltonian is best" is not a well-posed question here without
fixing the readout. — **ESTABLISHED**, n = 126.

### 1.2e THE EXACT LIMIT, with its scope stated

`tail_min ≡ argmin` in all 9 (H, α) cells, max difference **0**. That is an identity, not a
finding, and it is reported as one. Its content is this: Q6 (EXACT) puts every CVaR-optimal law on
the α-tail of `H`; the α-tail of `H` contains `argmin_P H`; the deployed readouts are order-based.

> **A pool-restricted CVaR-VQE whose training objective and readout energy are the same `H`, read
> out by an order-based rule, cannot select anything `argmin` cannot.** — **EXACT**

**Scope, because the claim is worthless without it.** It requires (i) pool restriction, (ii) the
same `H` for training and readout, (iii) an order-based readout. It **fails** if the readout
averages — and then the credit belongs to averaging, which is §1.2c — or if the readout energy
differs from the training energy.

**And it does not need convergence.** `OP.Field.best_z` is `argmin` over *everything the arm ever
evaluated*, so for any law, converged or not, a pool-restricted order-based readout returns
`min_{x ∈ S} H(x)` over the sampled multiset `S ⊆ P`, which is bounded below by `argmin_P H`. So
the bound holds for the whole trajectory, not just the fixed point, and CVaR training can only
change **which** members get sampled — never the attainable value. On a pool of K = 500 that bound
is reached by **500 evaluations of exhaustive enumeration**, and Sprint 20 measured that spending
**8192** evaluations does worse (Q3). That is Q3 with a proof attached, for the pool-restricted
case.

**Where the live question actually is, stated so this is not read as more than it is.**
`s20.qb2_opt.arm_vqe` — the selector Workstream A's mandatory matrix declares
(`"readout": "argmin of the training H over everything seen (Field.best_z)"`) — samples
**continuous torsions from per-residue basin mixtures**, not pool members. It is **not**
pool-restricted, so **D4 does not bound it.** The unbounded, live cases are exactly three:
(i) a generator that leaves the pool, (ii) a readout energy different from the training energy —
which A's `x_global` cross-readout cells measure — and (iii) an averaging readout, where the credit
belongs to averaging. `tail_member` in §1.2b prices the third end of the range: a uniformly random
member of the tail, i.e. the worst case for a converged law with a one-shot readout.

## 1.3 D1 — IS THE ENCODING LEVER A STEP-COUNT EFFECT?

### 1.3a Audit, before any rerun

**Two of the twelve cells are one measurement.** `AMB|nelder` and `AMBc|nelder` are
**bit-identical** — same `rmsd_ORACLE`, same `best_std`, same `iters` — on all 10 targets × 2
seeds, in **both** encodings, verified elementwise in `s20/results/qb2_enc.json` and
`qb2_opt.json`. `s20/qb2_lib.py:13–19` states the reason itself: `AMBc = sign(E)·log1p(|E|)` is a
strictly monotone transform, and Nelder–Mead is comparison-only. So is `best_of_N`, and so is any
order-based arm.

The identity was documented; the *count* was not corrected for it. Sign test over 12 cells with 11
negative gives p = 0.0063; over **11 distinct** cells with 10 negative it is **p = 0.0117**, and
**two**, not three, distinct cells have CIs excluding zero. — **EXACT** (the identity),
labelling correction (the count).

**Action for Workstream B, standing:** exclude every comparison-only arm from any AMB-vs-AMBc
conditioning table. Those rows are exactly zero by construction; they are identities, not nulls.

**The measurement itself survives the correct unit.** Recomputed at the brief's mandated unit —
mean Δ per target over the 11 distinct cells, paired bootstrap over 10 targets:
**−0.2167 [−0.3239, −0.1095], 9W/1L, sign p = 0.0215.** I do not dispute that something moved.

**The embedding takes fewer optimisation steps in every cell**, from the artefacts' own `iters`:

| arm | embedded | θ | ratio | ΔRMSD (AMB / AMBc / LEG) |
|---|---|---|---|---|
| `adam_fd` | 5.1 | 10.5 | **0.49** | −0.287 / −0.369 / −0.191 |
| `spsa` | 159.2 | 207.6 | **0.77** | −0.247 / −0.649 / −0.131 |
| `lbfgs_fd` | scipy, runs to its own convergence | | | −0.048 / −0.125 / **+0.073** |
| `nelder` | scipy (and the duplicate) | | | −0.272 / −0.272 / −0.138 |

Sprint 20 declared the halved-iteration confound for the FD arms and then wrote *"It does not touch
the RMSD column, which is simply the best structure found within the same 512 evaluations from the
same starts."* On this instrument that inference runs backwards: the start is retrieval-derived and
good, and optimising the deployed objectives harder makes the **structure** worse (Q3;
`search-saturates-discrimination-binds`). The embedding's own objective column agrees — it reaches
a **worse** objective in 9 of 12 cells. Supporting sign: `lbfgs_fd`, the arm least sensitive to
step count, carries the three smallest effects and the only positive one.

**And SPSA's deficit is a DIAGNOSTIC, not the declared dimension effect.** `arm_spsa`'s
`track_quality` spends three central-FD probes of `2·d` budget units; `d` doubles in the embedding,
so the probe consumes **288 of 512** units embedded against **144 of 512** in θ. A diagnostic
instrument was consuming 28% versus 14% of the two arms of the headline comparison.

### 1.3b The soundness gate

`arm_spsa_rec` / `arm_adam_rec` in `s21/d_enc.py` are verbatim copies of `s20.qb2_opt.arm_spsa` /
`arm_adam` with recording added and nothing else changed. Run with `track_quality=True`, the same
seeds, starts and budget, they must reproduce Sprint 20's **stored** `rmsd_ORACLE`, `best_std` and
`iters` bit-for-bit.

    GATE: 48 comparisons -- 3 targets x 2 kinds (LEG, AMBc) x 2 arms (spsa, adam_fd)
          x 2 seeds x 2 encodings -- 0 mismatches.  PASS.
    s21/results/d_enc_gate.json + .COMPLETE   (n_comparisons=48 mismatches=0)

The gate **fired 48 times** (a gate that never fires is not evidence).

**A defect of my own, found while tidying and fixed rather than explained.** The first gate run
printed **32/32, 0 mismatches** to the console and was then killed by `core.amber.memory_guard` at
95% physical memory **before it wrote its artefact** — so for several hours this document cited a
number with no file behind it. That is the vacuous-flag hazard wearing different clothes, and it is
exactly what I would flag in another lane. The gate was rerun on the **full default configuration**
and the figure above is read off the artefact, not off my memory of a console.

### 1.3c The grid, LEG half. n = 10, **4 seeds**, `track_quality` OFF, COMPLETE

**The mechanism is confirmed.** With the diagnostic off, SPSA's iteration counts are **exactly
equal** in the two encodings — 128/128, 256/256, 512/512 at B = 256/512/1024. Sprint 20's entire
159-vs-208 SPSA step deficit was the diagnostic. `adam_fd` keeps a 0.49–0.50 ratio, which is the
genuine 4n-vs-2n dimension effect and cannot be removed at fixed budget, which is why a budget grid
was measured instead of a single handicapped point.

**The gauge radius drifts, exactly as predicted.** `‖u‖` starts at 1.000:

| arm | B | `r_last` | `r_max` |
|---|---|---|---|
| `spsa` | 256 | **2.673** | 3.788 |
| `spsa` | 512 | **1.941** | 3.058 |
| `spsa` | 1024 | 1.747 | 2.827 |
| `adam_fd` | 1024 | 1.173 | 1.618 |

`arctan2` is scale-invariant, so the objective's gradient has **exactly zero** radial component and
every unit of radial motion is injected finite-difference noise. The effective angular step decays
as `1/‖u‖`, so the embedded SPSA arm is running an **unintended annealing schedule** worth roughly
a factor of two. — **ESTABLISHED**

**The outcome on LEG.** Negative = the embedding wins.

| comparison | `spsa` | `adam_fd` |
|---|---|---|
| A — matched **budget** (Sprint 20's comparison) | +0.215 [−0.055, +0.507] 5W/5L | −0.125 [−0.301, +0.038] 7W/3L |
| B1 — θ handicapped to the embedding's step count | +0.215 (steps already equal) | −0.052 [−0.139, +0.040] |
| B2 — embedding raised to θ's step count | +0.215 | −0.071 [−0.225, +0.071] |

Pooled at the target level over both arms: **A +0.045 [−0.109, +0.211], W/L 5/5, sign p = 1.00**;
B1 +0.082; B2 +0.072.

Sprint 20's LEG column was `adam −0.191`, `spsa −0.131`, both ns. **Mine flips SPSA's sign and
shrinks adam's to inside the MDE.** The LEG effect does not replicate. — **NOT MEASURED** on LEG.

**And the bar has to be the right one.** By §4.5, this comparison's own 80%-power MDE is
**0.242 Å**, not the brief's 0.084 — `sd(paired) = 0.273` at n = 10. So the honest statement is
that **my LEG null is uninformative below 0.242 Å**, which is weaker than the flattering version
and is the one that stands. It does not change the disposition (the point estimate is on the wrong
side for the lever, and Sprint 20's own LEG cells were ns), but it does mean this half **cannot**
rule out an encoding effect of a size the sprint would care about. Only the AMBc half can, and its
Sprint-20 counterpart has an MDE of 0.801 Å — see §4.5.

**A limitation I will not hide, and it cuts against my own hypothesis.** Panel C measures the
step-count response *within* each encoding (B = 1024 minus B = 256): `spsa/θ` −0.164 [−0.385,
+0.044], `adam/θ` +0.134 [−0.170, +0.526]. Both span zero and they point in **opposite**
directions. So on LEG I have established that **the confound exists** and that **the effect does
not replicate**; I have **not** established the direction of my proposed mechanism. That is
**NOT MEASURED**, and the AMBc half — where Sprint 20's two surviving significant cells live
(`AMB/spsa −0.247`, `AMBc/spsa −0.649`) — is the arm that decides it.

### 1.3d The decisive cell: AMBc / spsa. n = 10, **4 seeds**, COMPLETE

`s21/results/d_enc_AMBc_spsa_B512.json` + `.COMPLETE` (80 cells: 10 targets × 2 encodings ×
4 seeds, `track_quality` OFF, budget 512). This is the cell that matters — `AMBc/spsa` is the
largest Sprint-20 encoding effect (**−0.649 [−1.205, −0.165], "sig"**) and one of only **two
distinct** cells there with a CI excluding zero, the other having been the Nelder duplicate.

**The full budget grid was abandoned deliberately and the partial is preserved.** It was running
at ~25 min/target under four-way contention (52 ms per AMBER single point against a 6 ms nominal),
would have taken four hours, and was starving `tailprice` — which logged *"free 0.42 GB < 1.8;
waiting"* for long stretches. `_PARTIAL_d_enc_AMBc_fullgrid_3targets.json` is kept. The reduced
design is the *right* cut, not merely the affordable one: with `track_quality` OFF, SPSA's step
counts are **exactly equal** in the two encodings, so the matched-budget comparison **is** the
matched-step comparison and no grid is needed.

| | θ | emb |
|---|---|---|
| iterations | **256.0** | **256.0** — exactly equal |
| `disp_best` (rad, torsion space) | 5.012 | **4.492** — ratio **0.90** |
| gauge radius `‖u‖` | 1.000 by construction | **2.498**, max **4.027** |
| RMSD (built chain) | 5.107 | 4.812 |

    A  matched budget                       -0.295 [-0.636, +0.022]   W/L 6/4   med -0.194
    B1 matched steps (identical -- 256/256) -0.295 [-0.641, +0.027]   W/L 6/4
    B2 matched steps (identical -- 256/256) -0.295 [-0.631, +0.022]   W/L 6/4
    target-level pooled, sign p = 0.7539            own MDE ~ 0.53 A

> **Sprint 20's −0.649 [−1.205, −0.165] "sig" becomes −0.295 [−0.636, +0.022], a 55% smaller
> point estimate with a CI spanning zero, once the diagnostic's budget asymmetry is removed and
> the seed count is doubled.** My pre-registered falsifier — *the embedding still winning by more
> than the MDE with a CI excluding zero in both matching directions* — **does not fire.**

**What I will not claim.** Two mechanisms are jointly sufficient to explain the shrinkage — the
`track_quality` budget asymmetry (§1.3a) and Type-M magnitude exaggeration in an underpowered panel
(§5.1: post-hoc power **0.62**, exaggeration **1.27×** for that cell) — and **this design cannot
apportion between them.** And the panel's own MDE is **~0.53 Å**, so it **cannot exclude** an
effect of the size originally reported. The honest label for AMBc is **NOT MEASURED**, with a point
estimate 55% smaller than the one it replaces.

**What survives, and it is a real defect rather than a result.** At *exactly equal step counts* the
embedding still moves **10% less far** in torsion space (4.492 vs 5.012 rad) while `‖u‖` grows to
**2.498 (max 4.027)** along a direction whose true gradient is exactly zero. That is channel (b) of
the hypothesis — the **gauge-driven step-size decay** — still operating with channel (a) removed,
and it is consistent with the residual −0.295. **`_EmbField` has a bug, not a feature**: the
embedded arm injects finite-difference noise into a pure gauge direction and thereby anneals its
own step size by roughly a factor of two. Any future embedded run should project the update onto
the tangent space or renormalise `u` each iteration; until it does, an embedded arm and a θ arm are
running different optimisers, not different coordinates.

### 1.3e The powered replication: LEG and DIST at n = 60, **COMPLETE**

`s21/results/d_enc_LEG-DIST_n60_spsa-adam_fd_B256-512.json` + `.COMPLETE`: **60 targets × 2 kinds ×
2 arms × 2 encodings × 2 budgets × 4 seeds = 3840 cells**, `track_quality` OFF.

**Why this exists.** §5.1 shows an n = 10 panel's own 80%-power MDE on this comparison is
0.24–0.80 Å, so **neither Sprint 20's panel nor my first one could ever have resolved the encoding
question.** On the objectives cheap enough to afford it, the tuning instrument's first 60 targets
are used. This is **not paired to Sprint 20's panel**; it is an independent, better-powered
replication and is labelled as one.

Negative = the embedding wins. Fold-clustered CIs.

| cell | steps emb/θ | **A** — matched budget | **B1** — matched steps |
|---|---|---|---|
| LEG / `spsa` | 256 / 256 | +0.021 [−0.084, +0.128] 31W/29L | identical (steps equal) |
| **DIST / `spsa`** | 256 / 256 | **−0.005 [−0.038, +0.029]** 33W/27L | identical |
| LEG / `adam_fd` | 5 / 10 | **−0.130 [−0.204, −0.058]** 38W/22L | **−0.039 [−0.090, +0.010]** |
| DIST / `adam_fd` | 5 / 10 | **+0.091 [+0.024, +0.158]** 24W/36L | **+0.004 [−0.032, +0.039]** |

**Target-level pooled, the brief's mandated unit:** `A` **−0.0061 [−0.0448, +0.0329], W/L 30/30,
sign p = 1.0000**; `B1` −0.0049 [−0.0379, +0.0286]. CI half-width **0.039**, own MDE ≈ **0.055** —
**below the 0.084 bar. This is a POWERED flat null, not an underpowered one.** 30W/30L is as dead
as a comparison gets.

At matched **budget** the `adam_fd` effect is significant **in both directions** — −0.130 on LEG,
**+0.091 on DIST**. A change of coordinates that provably does not change the physics cannot help on
one energy and hurt on another. A change in **step count** can, and §1.3f shows it does, to 0.004 Å.

### 1.3f THE MECHANISM, quantitatively — and it refutes my own stated hypothesis

For `adam_fd`, the embedding at B = 512 takes the **same number of steps** as θ at B = 256. So the
objective's own step-count response, measured **in θ**, predicts the encoding effect as *minus
itself*. That is a prediction with no free parameter, testable on both objectives:

| | step response `θ(512−256)` | **predicted** effect | **observed** effect | **residual** |
|---|---|---|---|---|
| LEG / `adam_fd` | **+0.0912** [+0.078, +0.109] | −0.0912 | −0.1301 | −0.0389 [−0.122, +0.025] |
| DIST / `adam_fd` | **−0.0865** [−0.168, −0.018] | +0.0865 | +0.0905 | **+0.0040** [−0.023, +0.027] |

> **Both residuals span zero, and DIST's is +0.004 against a predicted +0.087.** The encoding
> "effect" **is** the objective's own step-count response evaluated at the step deficit — on two
> objectives, with **opposite signs**. — **ESTABLISHED**, n = 60

**And it refutes my own stated mechanism.** I proposed that the confound acts because *"fewer steps
helps here"*. **It does not.** The step response is **significantly signed in both directions**
depending on the objective: LEG **+0.0912 [+0.0777, +0.1089]** (more optimisation is *worse*) and
DIST **−0.0865 [−0.1678, −0.0179]** (more optimisation is *better*), **both CIs excluding zero**.
So *"optimising harder hurts"* is **not a property of this instrument** — it is a property of the
**Legacy** objective, and the **deployed distogram objective goes the other way**.

That is a correction to how I had been quoting `search-saturates-discrimination-binds`, and it
matches the brief's own §5 formulation, which I should have used from the start: *"the budget trap
is a property of BAD OBJECTIVES: searching harder hurts on a bad one, is neutral on a mediocre one,
and HELPS on a good one — always state that condition."* My data now measures that with CIs on
**both** signs. **O2 is CLOSED, and not in my favour**: the confound is fully explained, and the
mechanism is not the one I proposed.

### 1.3g Disposition of D1

**The encoding lever does not survive the step-count control, and the reason is now mechanistic
rather than merely statistical.** Every cell in which it was significant either (a) was an
identity — `AMB|nelder ≡ AMBc|nelder`, a comparison-only arm on a monotone transform (§1.3a); or
(b) had unequal step counts, loses significance when they are matched, and has its magnitude
**predicted to 0.004 Å** by the objective's own step response (§1.3f); or (c) shrinks 55% with a CI
spanning zero once a *diagnostic's* asymmetric budget is removed (§1.3d). Where step counts are
exactly equal the pooled effect is **−0.006 Å at 30W/30L with a powered MDE of 0.055**.
— **NOT SUPPORTED.**

Recommended disposition for `s21/BRIEF.md` §5 item 3: from *"determine whether this is robust and
exploit it"* to **NOT SUPPORTED**, and build no ansatz or representation on it.

**What I still have NOT shown.** That **no** encoding effect exists anywhere: the AMBc panel's own
MDE is ~0.53 Å and cannot exclude an effect of the size Sprint 20 reported. The AMBc label stays
**NOT MEASURED** with a 55%-smaller point estimate; the LEG+DIST label is a **powered null**. And
the durable engineering output remains the **gauge defect** of §1.3d (D16), which is a bug, is
fixable in one line, and will recur in any manifold-embedded arm this programme writes.

## 1.4 D7 — the brief's opening row, traced. My hypothesis REFUTED; a provenance correction stands

`s21/d_distobj.py` → `d_distobj.json` + `.COMPLETE`, n = 126, zero AMBER.

`s21/BRIEF.md` §2 opens the sprint with `argmin RMSD AMBER 4.990 | Legacy 5.487 | distogram 3.676,
pool MEAN 4.739` and builds the "argmin is not CVaR" framing on it. Traced to
`s20/agentB_FINDINGS.md` §6d, those four numbers are **n = 20 targets**, on the ideal-geometry
**rebuild** of each window's torsions, and the "distogram" column is `s20.qb2_lib.Ham("DIST")` =
`Σ_p (d_p − d̂_p)² / sd_p²`, a **squared-distance** functional — **not** the deployed selector,
which is `s12.instrument.shipped_score`, a **Bayes-risk** score off the full posterior. I
hypothesised the functional was the difference.

**It is not. My hypothesis is REFUTED, at n = 126 on one common basis:**

| | mean | median |
|---|---|---|
| `argmin` by the **shipped Bayes-risk** score | **3.4540** | 3.4779 |
| `argmin` by the **squared-distance** functional | 3.5040 | 3.4715 |
| **`ρ(E_bayes, E_sq)` per target** | **0.9570** | **0.9729** |
| `ρ(E_bayes, true)` / `ρ(E_sq, true)` | 0.5678 / 0.5621 | 0.6935 / 0.7003 |

`argmin_sq − argmin_bayes = +0.0500 [−0.0185, +0.1554]` fold-clustered, against that comparison's
own MDE of **0.1235**. My pre-registered falsifier fires: **the two distance functionals are
interchangeable**, which is good news — the matrix's `Distance` row is robust to which one is used.

**What does stand is a provenance correction, and it matters for the matrix.** On Sprint 20's own
`subset(20)`, window basis, I get `pool_mean 4.7318` and `pool_best 1.8882` against its reported
`4.739` and `1.901` — **agreement to 0.01 Å**, so the subset and pool reproduce and the
rebuild-vs-window basis shift is negligible (Workstream A measures it at **+0.011 Å** per member).

> **The 3.676-vs-3.454 gap is the SUBSET.** `subset(20)` is simply harder than the instrument:
> pool mean **4.732 against 4.453**, pool best 1.888 against 1.711. Anyone benchmarking an n = 126
> arm against the brief's 3.676 / 4.739 is comparing against a **harder 20-target subset**, and the
> shipped selector's own n = 126 argmin is **3.4540** — `s12/instrument`'s pinned constant.

**Ask:** annotate §2's row with `n = 20, subset(20), rebuild basis, squared-distance functional`,
or replace it with the n = 126 window-basis row `argmin disto 3.454, pool mean 4.453`. Either is
fine; the unattributed version is the fifth basis/provenance hazard this project has met.

## 1.5 D8 — auditing the coordinator's exhaustive-latent arm, which took up §3.1a's O7

`s21/latentfull.py` → `latentfull.json`, **75 targets** (every one with `2^n ≤ 8192`), all `2^n`
bitstrings decoded at their von Mises modes, **COMPLETE**. This is §3.1a's recommended control,
built and run by the coordinator within the hour. Audited here as a live claim.

**The primary carried an unstated functional mixture, and D7 is what found it.** The latent arms
were scored with `s19.qb_lib.Obj.raw` (the **squared-distance** functional) and the pool arm with
`I.shipped_score` (the **Bayes-risk** score) — exactly the mixture §1.3g measured at n = 126. The
artefact now carries both columns, so the matched primary can be computed from it directly:

| primary | mean | CI95 fold-clustered | W/L | own MDE |
|---|---|---|---|---|
| latent argmin (**squared**) − pool argmin (Bayes) | +0.2152 | [+0.0876, +0.3627] | 32/43 | 0.224 |
| **latent argmin (Bayes) − pool argmin (Bayes)** | **+0.1009** | **[−0.0234, +0.2343]** | 28/47 | 0.196 |

**Matching the functional halves the gap and its CI comes to span zero.** On the same enumerated
structures `Bayes − squared = −0.1143 [−0.2142, −0.0032]`, **median exactly 0.0000** — the two
agree on the typical target (§1.3g's median ρ = 0.973) and differ on a tail, and the tail pointed
in the direction that flattered the conclusion. Recommended label for the matched primary:
**NOT MEASURED**, +0.101 against its own MDE of 0.196, per-fold +0.13/+0.20/+0.33/−0.02/−0.09.
The coordinator's own framing — *"ties on the typical target, loses on a tail, wins nowhere"* —
survives the correction and is the right strength for it (28W/47L at a median of +0.000).

**The load-bearing numbers use one functional on both sides and are untouched.** Verified from the
artefact:

| | mean | CI95 fold-clustered | W/L | own MDE |
|---|---|---|---|---|
| `argmin − ORACLE` inside the latent (**squared**) | **+1.7127** | [+1.4056, +2.1887] | **0/71** | 0.436 |
| `argmin − ORACLE` inside the latent (**Bayes**) | **+1.5984** | [+1.2736, +2.0090] | **0/69** | 0.391 |
| `argmin −` the space's own **mean** | −1.1156 | [−1.4393, −0.7699] | 65/10 | 0.396 |
| latent ORACLE ceiling − the **incumbent** | −1.0822 | [−1.3077, −0.9120] | 62/13 | 0.395 |
| modes vs stochastic draws (the mode restriction) | −0.0830 | [−0.3260, +0.1564] | 43/32 | 0.255 |

Per-fold on the discrimination gap: **+1.60 / +2.56 / +1.63 / +1.54 / +1.30 — 5/5 same sign, and
0 of 71 targets go the other way**, under either functional.

> **The objective has real ordering skill — it beats the mean of the space it searches by 1.12 Å.
> The latent CONTAINS structures 1.08 Å better than the incumbent. And given EXHAUSTIVE access to
> every one of them, the objective still misses by 1.6–1.7 Å.** On the enumerable half of the
> instrument, `search-saturates-discrimination-binds` stops being an inference from six instruments
> and becomes a **property of the instrument**. The mode-restriction arm is what makes the
> exhaustiveness argument valid, and it costs nothing (−0.083, NOT MEASURED).

**Two defects in the module, one of which matters more now that the effect has halved.**

* **It writes its JSON non-atomically** — `json.dump(..., open(path, "w"))` with no `tmp` +
  `os.replace`, unlike `d_enc._save`, `tailprice`, `a_matrix._write` and `qb2_run.ck_save`, which
  all do it properly. **I caught a 45-row partial mid-write and computed a 5-target statistic on it
  before noticing** the file had held 75 rows a minute earlier. One line, and it silently produces
  wrong tables in a concurrent lane.
* **A second, unstated basis difference sits in the same primary.** The `latent_*` arms are **built
  chains** (`ob.raw` → `build_ca` from torsions) while `pool_argmin` is `W[order[0]]`, the retrieved
  window's **own coordinates** — window basis, not a rebuild. Workstream A prices the
  rebuild-minus-window shift at **+0.011 Å** per member, so the bias is small and points **against**
  the latent arm. At +0.215 that was 5% of the effect; at the matched **+0.101 it is 11%**. Fix:
  score the pool arm on `I.build_ca` of the window's own torsions — A's `_RB` arms already do this.

### 1.5a COMPLETED TO n = 126, and scored against predictions registered before looking

`s21/results/latentfull_hi.json`, **51/51, `complete: true`** — the `n = 14, 15, 16` targets
(16,384 / 32,768 / 65,536 configurations each, 2.13 M total). Merged with the 75 already done, the
enumeration now covers **the whole instrument**, and the length-defined subset caveat is gone.
Predictions P1–P3 are in `PREREG_D.md`'s D8 addendum, registered **before** this ran.

| merged, n = 126, fold-clustered | mean | CI95 | W/L | own MDE |
|---|---|---|---|---|
| discrimination gap `argmin − ORACLE` (**squared**) | **+1.8149** | [+1.6129, +2.0594] | **0/122** | 0.340 |
| discrimination gap `argmin − ORACLE` (**Bayes**) | **+1.7202** | [+1.4991, +1.9752] | **0/119** | 0.326 |
| `argmin −` the space's own **mean** | −1.3674 | [−1.5333, −1.2100] | 109/17 | 0.352 |
| latent **ORACLE** ceiling − the **incumbent** | −1.1766 | [−1.2935, −1.0645] | 104/22 | 0.325 |
| **P3 primary**: `bayes argmin − pool_argmin_rb` | +0.0851 | [−0.0176, +0.1691] | 53/73 | 0.197 |

> **Not one target in 126 has the objective's exhaustive argmin beating the latent's own ORACLE
> best.** 0/122 squared, 0/119 Bayes.

**P2 — CONFIRMED, and stronger than registered.** I predicted "0 or near-0 of 51 reversed"; it is
**0/51** on the new targets and **0/122** merged, per-fold signs uniform throughout.

**P3 — CONFIRMED, on the targets where it was least obvious.** The 51 are the **only** targets whose
latent exceeds the deployed 8192-evaluation budget — the only ones where the sampler was genuinely
*searching* rather than resampling (§3.1a). There the matched primary is **+0.1738 [−0.0079,
+0.3715], 19W/32L**: exhaustive search still does not beat the pool, and if anything does slightly
*worse* than on the enumerable-within-budget half (+0.0247). Cross-range difference +0.1491
[−0.1505, +0.4549] — NOT MEASURED. **Budget-versus-latent-size was not the operative variable.**

**P1 — NOT MEASURED, and my mechanism is REFUTED.** Two separate things.

* *The direction was right and inside my guessed band*: gap **+1.9651** on `n ≥ 14` against
  **+1.7127** on `n ≤ 13`, a difference of **+0.2524** (squared) / +0.3009 (Bayes), against my
  registered guess of +0.2 to +0.6. But the cross-range difference **CI spans zero**
  ([−0.2343, +0.7391]). **NOT MEASURED** — the label I pre-committed to whichever way it landed —
  and now for two independent reasons: the min-of-N confound I declared in advance, *plus*
  insufficient power on the difference itself.
* *My stated mechanism is contradicted by the data.* I predicted `latent_oracle` would **improve**
  with `n` by min-of-N, its set growing 128× from `n = 9` to `n = 16`. It goes the **other way**:
  `corr(n, latent_oracle) = +0.297`, per-length means rising **1.475 → 2.025**. The harder-target
  effect swamps min-of-N entirely. The arm that actually moves with length is the space's **mean**:
  `corr(n, latent_mean) = +0.677`, against +0.185 for `argmin`, +0.297 for ORACLE and +0.144 for
  `pool_argmin_rb`. **Right direction, wrong reason, and unmeasurable either way.**

That is the **third** of my hypotheses this sprint whose *mechanism* was wrong even where its
direction survived — after D1's *"fewer steps helps here"* (§1.3f) and R1's *"diversity drives
averaging"* (§4.1). The pattern is worth naming: **a directional prediction can survive while the
story attached to it fails, and only a mechanism-specific measurement separates them.** In all three
cases the separating measurement was cheap and the mechanism was the part that was wrong.

**A limitation I raised on my own task and the coordinator fixed rather than labelled.** `min`-over-
`2^n` also means `latent_oracle` at large `n` is partly measuring **set size**, which cuts against
the *"the latent contains the answer"* reading at large `n`, not only against the gap. My
`corr(n, latent_oracle) = +0.297` shows the ORACLE is not artificially good in **absolute** terms —
but that is a *level* argument, not a decomposition, and it is **not** a substitute for the
matched-set control. `s21/latentrank.py` now carries an `oracle_8192` column: the minimum over a
fixed, seeded 8,192-configuration subset, identical in size at every `n`, coinciding with
`latent_oracle` by construction at `n ≤ 13` — which is both the control and its own wiring check.
Until that lands, **the cross-`n` magnitude stays NOT MEASURED.**

### 1.5b THE TOP-M LADDER WAS min-of-M, AND THE HAZARD I REGISTERED OVERTURNED THE HEADLINE

`s21/latentrank.py` reported an ORACLE ceiling of **1.986 Å inside the objective's top-512**
against the deployed argmin's 3.281 — read as *"a reranker over the top-512 reaches below the
mission target."* Before its numbers existed I registered (`PREREG_D.md` D9) that **the top-M
ceiling is a minimum over M**, so a curve that falls from M = 1 to M = 512 falls *partly because
512 draws beat 1 draw*, and that the matched control is a **min-of-M over a RANDOM window of the
same M** — the vertical gap between the two curves being the ordering skill in the currency a
reranker spends. The coordinator implemented it. Measured, n = 75, `GAP = rand_M − top_M`
(negative = the objective helps):

| M | top-M | rand-M | **GAP** | CI95 | W/L |
|---|---|---|---|---|---|
| 1 | 3.281 | 4.466 | **−1.184** | [−1.448, −0.918] | 63/12 |
| 8 | 3.068 | 3.270 | −0.202 | [−0.490, +0.074] | 34/41 |
| 64 | 2.498 | 2.459 | +0.039 | [−0.204, +0.281] | 39/36 |
| **512** | 1.986 | **1.879** | **+0.107** | [−0.051, +0.300] | 47/28 |

> **The objective's ordering skill is concentrated entirely at M = 1 and is exhausted by M ≈ 8.**
> At M = 1 it is enormous — the argmin is 1.184 Å better than a random configuration, 63W/12L. By
> M = 64 it is zero. At M = 512 the objective-ordered set is if anything **worse** than 512 random
> draws, which reach **1.879** and beat the incumbent by −0.886 [−1.175, −0.614] where the ordered
> 512 manages −0.779. **The entire fall of the ladder from 3.281 to 1.986 is min-of-M.**

**And my own registered prediction on it is REFUTED.** I wrote that the gap would be *"positive and
substantial at every M"*. It collapses to zero by M = 64 and **reverses sign** at M = 512 — I
predicted the objective ordered better than it does, which is the comfortable direction. Scored as
**REFUTED**. What survives is only that the **control was necessary**, which the numbers establish
better than the prediction could have.

**What the numbers support, which is less comfortable than the ladder read alone.** *"Dig deeper
into the objective's ranking"* is **dead** — there is nothing down there its ordering knows about.
Sharper than "the skill is exhausted": **the objective is a good FILTER and a useless RANKER.** It
knows which single configuration is best far better than chance and knows nothing about the
ordering of the next 500. Any architecture spending budget on the *ranking* is spending it on the
part the objective does not have. What is alive is that **the latent is rich enough that 512
arbitrary draws contain a sub-2.0 Å structure** — and any route to that richness is a **different
selector**, not a deeper readout of this one.

**Reconciling this with Q1's 5.62th percentile,** because they look contradictory. The ORACLE-best
sits at median percentile **5.62** — genuinely better than the 50.0 null, so the ordering does carry
information — but the mean is **20.27** with a max of 99.91, and a top-512 filter on a 4096-config
latent keeps 12.5%. **Median rank skill and top-M containment are different quantities, and only
containment is what a reranker spends.** Containment is governed by the **upper tail** of the
percentile distribution, not its centre. My Q1 was about the first and landed at the boundary
(registered 0.1st–5th, measured 5.62 — **outside my interval; I score it a miss, and only the
reasoning survives**). D9 was about the second, and it is what decides the architecture.

**My min-of-N confound on the ORACLE arm is now measured too, and it is REAL BUT TINY**: at
`n ≥ 14`, the ORACLE over a matched 8,192-configuration subset is **1.909** against **1.872** over
the full latent — **+0.037**. Correct in kind, negligible in size. Worth having measured rather than
assumed in either direction.

### 1.5f D10 — CAN A NATIVE-FREE READOUT TRAVERSE THE OBJECTIVE'S TOP-M ON THE LATENT?

`s21/d_lrank.py` → `d_lrank.json`. 75 targets (`n ≤ 13`, `2^n ≤ 8192`), every latent mode
enumerated, four readouts inside the objective's top-M against a **matched-count random M-window**,
`M ∈ {1, 8, 64, 512}`, `R_DRAW = 8`. Pre-registered as D10 with T1/T2/T3 and their falsifier;
the containment column and the best-of-K null were registered as D11 before either was read.
**FINAL TABLE AND DISPOSITIONS: see the artefact and the closing message — the run's own flag
requires all 75 rows to carry every `nf_*` and `ORACLE_contained*` key.**

*Interim at n = 45, labelled as interim, is what the coordinator was sent to prevent a withdrawn
prediction reaching the report:* **T1 SUPPORTED** (the averaging readout beats the deployed argmin —
`top64_avg` 3.163, `top512_avg` 3.265 against 3.449); **T2 REFUTED on both halves** (the
avg-vs-random gain is **−0.478** at M = 64 and **−0.327** at M = 512, against the pool's **−0.377**,
clearing its own MDE at every M — larger, not smaller, and my error-coherence mechanism with it);
**T3 SUPPORTED** (`top512_avg` 3.265 against `ORACLE_top512` 2.211, a 1.05 Å gap unclosed).

> **The row I did not predict, and the one worth keeping.** `top{M}_argmin` is 3.449 at every M by
> construction — the argmin of the top-M *is* the global argmin. `rand{M}_argmin` falls
> **4.665 → 4.040 → 3.684 → 3.494**, ending at **−0.045 [−0.177, +0.137], NOT MEASURED.**
> **The objective's argmin over 512 random latent draws is statistically indistinguishable from its
> argmin over the entire 8192-configuration latent.** Exhaustive enumeration buys essentially
> nothing over 512 draws for the deployed readout — a **third** independent route to search
> saturation, and the most operationally direct: not *"searching harder does not help the
> structure"* but *"the search is DONE at 512 draws"*, measured against its own exact optimum. It
> also prices the budget ladder: **the difference between 512 and 8192 evaluations is under the
> MDE.**

**And an independent cross-check of §1.5b in a different instrument:** `rand512_ORACLE` **1.975**
against `ORACLE_top512` **2.211** — the random window's ceiling is *better* than the ordered
window's, reproducing D9's reversal from a different construction on a different subset.

### 1.5c AUDITING `latentsel.py` — the fork verdicts I was asked for, computed rather than requested

The coordinator sent the **operator-fork list before reading the numbers** — the first application
of the rule in §2 — and asked my view on one undecided fork. All three rows I recommended were
derivable from their artefact, so I computed them rather than asking for a rerun. **Final, n = 126,
`complete: true`.**

| | mean | CI95 fold-clustered | W/L | own MDE |
|---|---|---|---|---|
| **PRIMARY** `lat_avg75 − pool_avg75` (basis-matched) | **+0.3898** | [+0.3225, +0.4793] | 47/79 | 0.242 |
| footnote `lat_avg75 − ship_avg75` (composite) | +0.4053 | [+0.3440, +0.4950] | 47/79 | 0.240 |
| **the basis price** `pool_avg75(rebuilt) − ship_avg75` | **+0.0155** | [+0.0093, +0.0232] | 46/80 | 0.029 |

(`ship_avg75` = **3.0483** — the incumbent point cloud, bit-exact for the third time in this
document, which is the soundness check on the whole comparison.)

**The verdict on the fork: basis-matched is the primary, but measure the basis price.** A *source*
comparison asks whether swapping the candidate source changes the answer with the operator held
fixed, so everything except the source must be identical, basis included. But `vs ship_avg75` is a
**composite** of source swap and basis change, and leaving it composite asks the reader to trust
which dominates. Measured, **the composite is 96.2% source and 3.8% basis** — and it needed
measuring rather than borrowing Workstream A's per-member +0.011, because the rebuild can change
**which member is the medoid**, and the medoid sets the frame the superposition is computed in.

**Hazard (b) was a different arm, not a fork, and the sensitivity closes it in a line.** On the 9
targets with `n = 9`, `2^9 = 512`, so "512 draws without replacement" is *deterministic* — the
exhaustive latent with zero variance, handed the ORACLE-richest possible set while every other
target sees 6–50% of its own. The label is wrong on 12% of the panel and wrong in the direction
that favours the latent. Pre-declared sensitivity:

    PRIMARY on all 126        +0.3898 [+0.3225, +0.4793]
    PRIMARY on n >= 10 (117)  +0.4266 [+0.3473, +0.5027]
    PRIMARY on n >= 12  (93)  +0.4745 [+0.3412, +0.6527]

Dropping them does not move the primary; **it enlarges it monotonically as the length floor rises**.
That is the D9 min-of-M lesson showing up in a third place: at `n ≥ 12` a 512-draw window is
≤ 12.5% of the latent, and the latent arm does worse the *less* of its own space it sees — while the
pool arm always sees 500 of ~500. **The latent's apparent competitiveness on short peptides is a
coverage artefact, not a property of the source.**

**What the primary says, stated plainly because negatives attract soft language: the source swap
loses.** The latent through the *identical shipped operator* is **+0.41 Å worse** than the pool, CI
excluding zero, effect above its own MDE, 44W/76L. And it is not only extraction: `lat_oracle`
2.278 against `pool_oracle` 1.718 says the latent's 512 draws **contain** worse structures than the
pool's 500 windows on this panel. The score halves the deficit (+0.86 score-blind → +0.41
score-filtered) and does not close it.

### 1.5d MY CONSISTENCY PREDICTION IS REFUTED, AND THE RECONCILIATION IS THE FINDING

I predicted `lat_avg75` and `lat_avg75_rand` would be **indistinguishable**, on the reasoning that
D9 had shown the objective's ordering worthless past its top handful — and I said that if they
differed materially, one of the two measurements must be wrong.

    lat_avg75 - lat_avg75_rand   -0.4290 [-0.5077, -0.3317]  91W/35L  MDE 0.191   (n = 126)

They differ by **0.43 Å**, 91W/35L at n = 126. **Neither measurement is wrong — my inference was.** D9 and `latentsel`
measure **different functionals of the same set**:

| | what it compares | result |
|---|---|---|
| **D9** | the set's **BEST** (min over M) | no skill past M ≈ 8 |
| **latentsel** | the set's **MEAN** (through the averaging readout) | **0.43 Å** of skill at M = 75 |

A score-filtered set can have a much better **mean** while having no better **minimum**, and the
terminal operator consumes the mean: `operator-consumes-set-mean` has it exactly —
`d_out = 1.16·d_set_mean + 0.04·d_set_best, R² 0.89`.

> **The objective's ordering is USELESS for finding the best member of a set and worth 0.43 Å for
> shaping the set's MEAN. Which of those matters is decided by the READOUT.** D9 and `latentsel`
> disagree only because they read the set out differently.

I therefore **withdraw** my earlier phrasing *"a good filter and a useless ranker"* in that
unqualified form. It is my own §1.2d — a Hamiltonian's disposition flips with the readout — arriving
for the **third** time, now inside the latent rather than the pool: **no statement about an
objective's usefulness is well-posed until the readout is fixed.**

### 1.5e THE FAILURE CELL'S SIGN IS FORCED BY THE SPLIT

The coordinator's L19 resolves the containment paradox: the objective contains the ORACLE-best in
its top-512 on **59%** of targets and is **worse than random** on the other 41%, and the aggregate
`+0.107` is the average of `+0.396` and `−0.851`. The bimodality is real and it is the right
mechanism. They correctly note the *contained* cell is a restatement of the split. **The failure
cell is not clean either, and it is the one called load-bearing.**

    C := 1{ORACLE-best is inside the objective's top-512}.
    Conditioning on C = 0 SELECTS targets where the top-512 provably misses the best config.
    It imposes NO corresponding condition on rand-512.

So the C = 0 comparison is between one arm **selected to have failed** and one arm that was not.
**The sign is forced** — any selector, including a perfect one, looks bad on the subset defined by
its own failures.

**What is not forced is the magnitude, and that is where the finding lives.** On the failure cell
the top-512 sits **1.23 Å** above the latent ORACLE (3.533 vs 2.305) while rand-512 sits only
**0.38 Å** above it (2.682). The conditioning forces `top-512 > ORACLE`; nothing forces it to be
three times further away than an arbitrary window. **That** is the concentrated-wrong-region claim
and it survives.

**The clean version needs no conditioning at all.** Report unconditionally over all targets
`ORACLE_top512 − latent_oracle` and `rand512_ORACLE − latent_oracle`; their difference *is* the D9
gap in excess-over-ORACLE terms, with no split and no selection effect. Use the bimodal
decomposition as **mechanism**, not as evidence, with the failure cell's sign marked as forced. If
the conditioned comparison is to carry weight, condition **both** arms on their own containment
failure and compare like with like.

**And this is the third arrival of one pattern this sprint** — D13 (the MDE is a property of the
comparison, not the instrument), D14 (a marginal correlation is a correlate, not a contribution),
and L19 (an aggregate null is the average of two opposite regimes). All three are *an aggregate
hides a structure that decides the disposition*, and all three were found by asking what operator
produced the number.

**Why §1.4 was worth an hour on a null.** Its useful output was never *"the functionals differ"* —
that hypothesis was **REFUTED**. It was that measuring whether they differ, on a common basis at
n = 126, is what let a **known-direction 0.1 Å bias be identified inside another lane's primary an
hour later**. A null with a measured magnitude is a calibration; a null without one is nothing.

## 1.6 D3 — WHICH MATCHED-DISPLACEMENT NULL? The "72%" is a ratio against a skilled control

`s21/d_null.py` → `d_null.json` + `.COMPLETE`. No structure rebuilt, no energy evaluated: every
number recomputed from Workstream C's own Sprint-20 artefact, target-paired, **fold-clustered CI
quoted with the i.i.d. one beside it**.

`s21/BRIEF.md` and my own commission both quote L12 — *"72% of AMBER's damage is move size, not
direction"* — and instruct me to hold the λ-continuation to "the matched-displacement null, the one
that matters". **There are two, they are matched to the same per-coordinate RMS torus magnitude
from the same starts, and they disagree about the sign of AMBER's direction.**

| AMBER, magnitude 0.577 rad/coordinate | mean | CI95 fold-clustered | W/L |
|---|---|---|---|
| minimiser − start | +0.6198 | [+0.3213, +1.0037] | 4/26 |
| **isotropic** null − start | **+0.9977** | [+0.7826, +1.2252] | 3/27 |
| **toward_member** null − start | **+0.4438** | [+0.3884, +0.4891] | 6/24 |
| minimiser − isotropic | **−0.3779** | [−0.7151, −0.0763] | 23/7 |
| minimiser − toward_member | +0.1761 | [−0.1248, +0.5741] | 12/18 |
| **the two nulls differ** | **−0.5540** | **[−0.8310, −0.2879]** | **24/6** |

    share of the move's cost reproduced by the null:   toward_member 72%   isotropic 161%

**The number is right and the label is wrong.** `s20/c_land_null.py:24`, `agentC_FINDINGS.md` §4
and `LEDGER.md` L12 each call `toward_member` *"realisable, native-free, **zero information**"*. It
is realisable and native-free. It is **not** zero-information: a fixed-size move toward a uniformly
random pool member is, in expectation, a move toward the **pool centroid**, and pool consensus is
the one native-free signal this programme has established has positive in-band skill (S12's
score-filter + consensus medoid, −0.172 [−0.316, −0.027]; memory
`consensus-is-the-only-in-band-discriminator`, `consensus-is-outlier-avoidance`). Here it is worth
**0.554 Å** against isotropic at AMBER's magnitude and **0.081 Å** at Legacy's, both CIs excluding
zero.

`s21/BRIEF.md` §7 rule 4 asks for a *plausible-but-uninformative* control rather than a uniform
one, and `toward_member` is exactly the right control. But "uninformative" and "zero information"
are different bars, and the difference is what "72%" means: the same arithmetic against the
isotropic null reads **161%**, i.e. AMBER's direction is **better** than noise.

**Workstream C reported both.** `agentC_FINDINGS.md` §4 contains the −0.3779 row and the sentence
*"AMBER is not moving randomly — it beats an isotropic move of the same size by −0.378. It is moving
too far."* The record is complete; the **headline** is the half that survived into the brief. This
is a labelling correction to a claim the originating lane stated correctly.

**Binding consequence for this sprint's λ-continuation, registered before C's artefact exists:**
every continuation arm must be reported against **both** nulls, and **beating `toward_member` is
the bar** — beating `rand_iso` is not, because AMBER already does that. And a continuation whose
gain evaporates against `toward_member` has lost to *consensus*, which is a different and more
interesting statement than losing to noise. — **ESTABLISHED**

---

# 2. PRIORITY 2 — THE AUDIT

**Benchmark seal.** Verified by hash, byte-identical to Sprint 20's certificate (top of file). Not
opened, not read, not derived from. Future sprints can check byte-identity against that line
without ever learning a target name.

**The mission ladder reproduces exactly, and its labels are now right.** From
`bench_results/cache/1fc9f2dcf489e2fb`, n = 126:

    rmsd_avg   3.0483  (median 2.8373)   POINT CLOUD    -- the brief's 3.048
    rmsd_fit   3.2041  (median 2.9661)   BUILT CHAIN    -- the brief's 3.204, the incumbent
    rmsd_arm   3.2148  (median 2.9661)   BUILT CHAIN
    rmsd_full  3.2355  (median 2.9757)   REPAIRED EMISSION -- the brief's 3.236

`s21/BRIEF.md` §1 labels all three correctly (point cloud / built structure / deployed emission).
**My predecessor's Sprint-20 correction stuck**, and the label drift that ran from Sprint 12 to
Sprint 19 is not present in this sprint's brief. Recorded as a negative audit result, which is what
it is.

**Completion flags.** One live defect found and reported (§1.1b). Every artefact I wrote requires
the full configuration and **removes** its flag when the condition fails:
`d_cvarop.COMPLETE` requires n = 126 with all 3 Hamiltonians × 3 α × 4 readouts × 3 controls finite;
`d_sim.COMPLETE` requires every register size the 126 targets use, at every tested depth, with the
normalisation check passing; `d_enc.COMPLETE` requires every target × kind × arm × encoding ×
budget × seed with a finite RMSD.

**CI construction, quoted with the CI.** `s12.instrument.paired` is i.i.d. over targets, not
fold-clustered, across 293 call sites (s20 X4). Every disposition-bearing interval in this document
carries **both**. In `d_null.py` I used `s18.phys_lib.paired`'s `ci_fold` as primary — Workstream C
had already caught this themselves in `s20/c_land_null.py:64–77`, which is the correct disposition
and is noted here so it is not re-discovered a fourth time.

**Gates and their firing counts.** `d_enc`'s soundness gate fired **32 times, 0 mismatches**.
`d_sim`'s normalisation check fired **64 times, all passing** (Σp = 1.0000000000). `d_cvarop`'s
`tail_min ≡ argmin` identity check fired **9 times**, max difference 0 — and it is reported as an
identity, not as a finding.

**Shared referent floor.** No correlation between two deviations-from-a-common-reference is claimed
in this document, so the 0.505 floor does not bind on anything here. Stated explicitly because the
absence of the check is otherwise indistinguishable from having forgotten it.

**Basis.** Every structural row above states point cloud / built chain. The one place it bites this
sprint is §1.2b: `tail_avg` is a **point cloud** and the other three readouts are real windows.
Reading the `tail_avg` row against the `argmin` row is a basis mismatch — the fifth instance in this
project — and the table carries a basis column so it cannot be made silently.

**Seeding.** `s15.seed.stable_rng` throughout. `d_enc` uses **4 seeds** on every variational arm
(§6 minimum; Sprint 20's encoding panel used 2, against a recorded 0.200 Å within-target seed sd
that is 2.4× the MDE — an n = 10 × 2-seed panel is thin for a 0.2 Å effect, and that is a second,
independent reason the Sprint-20 encoding cells were fragile).

**A defect of my own, in the same class but not the same failure mode.** I restarted `d_lrank`
twice to add arms, and each restart launched a new process while an older one was still alive —
**three processes writing one artefact**, each clobbering the others. I noticed only because the row
count went *backwards* between two reads (70 → 30).

> **Atomicity protects readers from a half-written file. It does nothing against a second writer.**
> Two atomic writers to one path produce a perfectly well-formed file containing whichever wrote
> last, silently mixing two configurations — here, one build *with* the containment columns and two
> *without*.

The fix is a lock or a **configuration-derived output path**, not stronger writes.
`s21/d_enc.py` gets this right, and by accident: its tag encodes kinds, arms and budgets
(`d_enc_LEG-DIST_n60_spsa-adam_fd_B256-512.json`), so two differently-configured runs cannot collide.
`d_lrank` did not. **What actually caught it was the completion flag demanding the full key set
rather than just the row count** — a mixed file could not have passed — which is the argument for
writing flags that way even when the row count seems sufficient. Recorded beside the non-atomic-write
note because the two look identical and are not: one is a **read** hazard, the other a **write
collision**, and the standard mitigation for the first does not touch the second.

**Two smaller live-artefact notes, reported because a flag and a hash are claims.**
`s21/results/_SMOKE_a_matrix.json` carries `complete: true` at `n_rows = 2`. The `_SMOKE_` prefix
does the real work and the row-level check in `a_matrix._required` is the best one in the sprint —
it demands every Hamiltonian × arm × seed × α × readout cell a row promises, which is exactly the
standard §6 asks for. But `complete` still compares against the `n` the run was *called* with, and
the same file records `N_TARGETS: None` in the config that `cfg_hash` is computed over — so the
hash does not pin the run's target count. One line each: require `n_expected == N_TARGETS` for the
flag, and record the `n` actually used in `CFG` before hashing. My predecessor self-caught the
identical smoke-flag hazard in Sprint 20; it is the recurring one.

**A hazard, not a live defect.** `s16.energy_lib.legacy_components_of_windows(seq, PHI, PSI)` and
`s13.qarch_lib.amber_energies(space, S)` take arguments in different orders and different spaces
(continuous torsions vs discrete lattice indices) while both being reached for as "the Hamiltonian
on the pool". One sprint has already lost an hour to it (§1.1a). A one-line assertion on the first
argument's type in each would have turned a silent 3-target skip into an immediate error.

---

# 3. PRIORITY 3 — LITERATURE, AND THE REACHABILITY VERDICT

## 3.1 The specific question: is a quantum-resource claim reachable here at all?

**Answered by measurement, not by citation.** `s21/d_sim.py` → `d_sim.json` + `.COMPLETE`.

The question as posed — *at what bond dimension does classical simulability break for this circuit
class* — has an answer that does not involve the bond dimension. From source,
`s20/qb2_opt.py:322`, the deployed circuit is

    Q.MPSAnsatz(n, layers=2, final_ry=True, entangler="cnot")    with  n = F.n = RESIDUES

so the register is **one qubit per residue** and the tuning instrument is 9 ≤ n ≤ 16. The whole
Hilbert space has 512 to 65,536 basis states. I enumerated the circuit's **complete** probability
distribution exactly, for every register size the 126 targets use, at layers 1, 2, 3, 4, 6, 8, 12
and 20 — up to and past the point where `χ = 2^L` saturates `χ_max = 2^⌊n/2⌋`, i.e. at **maximal
entanglement for this topology, ring closure included**:

| n | states | layers | χ | χ_max | route | time | Σp |
|---|---|---|---|---|---|---|---|
| 16 | 65,536 | **2 (deployed)** | 4 | 256 | MPS | 654 ms | 1.0000000000 |
| 16 | 65,536 | 8 | **256 (saturated)** | 256 | dense | 27 ms | 1.0000000000 |
| 16 | 65,536 | **20** | **256 (saturated)** | 256 | dense | **64 ms** | 1.0000000000 |
| all 8 register sizes × all 8 depths | | | | | | **17.7 s total** | all normalised |

The two routes are the point. The **MPS** route costs `O(n·χ⁴)` and therefore blows up with depth —
that is a property of the *representation*. The **dense statevector** route costs `2^n` and is
**independent of depth**.

> **(a) There is no bond dimension at which classical simulability breaks for this circuit.** Not at
> χ = 4, not at χ = 256, not at any depth. The register is too small. Sprint 20's Q7 (bond
> dimension 4, ≤16-state HMM) is **true and is not the binding constraint**. — **ESTABLISHED by
> measurement**

### 3.1a And the budget already swallows the whole space on most targets

`s19/qb_lib.py:198` — `draw_from_basins(bits, …)` takes `bits` of shape `(B, n)` with `n` =
**residues** and `bits ∈ {0,1}`: one qubit per residue, **two basins each**. So the entire latent
space is `2^n` configurations. Against the budgets this programme spends:

| | |
|---|---|
| \|latent\| over the 126 targets | min **512**, **median 8192**, max **65,536** |
| budget 512 | `2^n ≤ budget` on **9/126** targets (7%) |
| budget 4096 | on **52/126** (41%) |
| **budget 8192** (Sprint 20's) | **on 75/126 (60%)** |

> **At Sprint 20's 8192-evaluation budget the budget equals or exceeds the ENTIRE latent space on
> 75 of 126 targets, and the median target's latent space is exactly 8192.** On a majority of
> targets the deployed VQE is not searching a space it cannot enumerate — it is *resampling* one it
> could have enumerated within budget.

That changes what Q3 means. *"No sampler, quantum or classical, at 8192 evaluations beats the
zero-evaluation retrieval pool"* is **not** a statement that the search was too hard on those
targets: exhaustive enumeration of the latent was affordable and would not have helped, because
what binds is the **objective's ordering over the latent**, not the sampler's reach.
`search-saturates-discrimination-binds` now has an **exhaustiveness argument** attached, not only
six empirical instruments. — **ESTABLISHED** (arithmetic + source)

**Two things I would act on.** (i) Any budget-scaling arm on this latent must print `|latent| =
2^n` beside the budget; a sweep from 512 to 8192 crosses the point where the budget swallows the
whole space on 60% of targets, and the curve's shape past that point is a resampling artefact, not
a search curve. (ii) The cheapest available control for **any** sampler on this latent is
**exhaustive enumeration on the `n ≤ 13` targets** — 75 of 126, costing `2^n ≤ 8192` evaluations,
no more than the samplers already get. It yields the **exact argmin of the objective over the
latent**, which is the true optimum of the search half. If that does not beat the pool either, the
search half is closed **by construction** on those 75 targets rather than by six instruments.

Sharper still on parameter count: the deployed ansatz has `3n` = **27 to 48** parameters over
`2^n` = 512 to 65,536 outcomes. It is a **restricted parameterisation of a small classical
categorical distribution** — one that has `2^n − 1` free parameters and can be written down and
optimised directly. The circuit is not adding representational power; it is subtracting it.

**(b) Trainability of anything that escapes.** To reach a register where classical simulation is
genuinely hard you need n ≳ 50 qubits, i.e. ~50-residue chains, which is outside this benchmark
entirely. There the standard expressibility–trainability tension applies, and so does
**Cerezo et al., "Does provable absence of barren plateaus imply classical simulability?"**
(arXiv:2312.09121; Nature Communications, August 2025), which collects case-by-case evidence that
the structure removing barren plateaus is the same structure admitting classical simulation —
barren plateaus are a curse of dimensionality, and current fixes encode the problem into small,
classically simulable subspaces; a *soft dequantization*. Its stated caveats are real and I do not
overstate them: average-case arguments, smart initialisations, models outside their assumptions,
possible provable superpolynomial advantages, and the fact that a quantum device may still be
needed for an initial **data-acquisition** phase. **None of the caveats helps here**, because here
the data-acquisition phase is unnecessary: the law is tabulable in 64 ms.

**(c) A demonstrated separation on a comparable continuous-variable encoding: none found.**

> **VERDICT.** A quantum-resource or advantage claim is **NOT REACHABLE** on this instrument, for
> **any** ansatz, while the register is one qubit per residue on 9–16-mers. This does **not** touch
> Q4 — the VQE genuinely trains, and a 48-parameter restricted model can be a perfectly good
> optimiser. It says that nothing measured here can ever be *evidence of quantum resource*, and
> that ansatz redesign cannot change that. **If any part of Workstream B's night is a
> quantum-resource claim, that part is closed before it starts.**

## 3.2 The 2026 continuous-space preprint, read for its methods

**arXiv:2609.02113**, *"Logarithmic-scale variational quantum eigensolver for off-lattice protein
structure prediction in continuous torsional angle space"*, 2 September 2026. Three things, and
they are not the ones the abstract leads with.

**Its own numbers carry our basis problem, in one sentence.** *"Chignolin reached a 0.623 Å Cα RMSD
in retained snapshots and 1.199 Å in final models; Trp-cage achieved a 2.501 Å RMSD among snapshots
(3.512 Å in final models)."* That is a 1.9× and a 1.4× gap between a best-of-set and what the
pipeline emits — **the same object as this programme's ~1.9 Å selection gap**, reported as two
numbers in one sentence with no statement of how the retained snapshot was chosen native-free.

**It reproduces this programme's central negative, independently.** *"The custom energy function
performed best overall, though **energy-ranking imbalances persisted across sampled landscapes for
all functions**."* That is `nothing ranks within the pool` / `the objective does not rank the
native`, found by an external group, on different peptides, with different energies (Rosetta,
OpenMM and a custom hybrid). Counting Sprint 20's two, this is a **third** independent external
corroboration and it is the strongest of the three because it is stated as a limitation rather than
as a finding.

**The "O(log₂ N) qubits" headline is a trade, and the paper's own last sentence says so:** *"By
converting physical qubit constraints into circuit **depth** constraints."* Reading N torsions off
`2^q = N` basis-state probabilities does not compress information — it moves the cost into
measurement shots, and resolving N probabilities each of size ~1/N costs `O(N²/ε²)` shots.
Separately, the simulator readout — *"extracts molecular torsions from relative **phases** in
statevector simulations"* — is not measurable on hardware at all without tomography, which is why
the hardware runs use a different **CDF decoder**. Those are two different algorithms reported
under one name; the hardware figure (**best 1.758 Å** on `ibm_cleveland` / `ibm_miami`) is the only
one that prices the hardware path. **No control against an untrained circuit or a best-of-N sampler
is reported.**

## 3.3 QuPepFold, and a statement about the field that is worth more than either paper

*QuPepFold: a Python package for hybrid quantum-classical protein folding simulations with
CVaR-optimised VQE* (PLOS One, PMC12893577) is a fourth group running the same CVaR-VQE
construction on peptides. Read for its methods:

| | QuPepFold |
|---|---|
| model | **tetrahedral LATTICE**, coarse-grained, `H = H_geom + H_chir + H_steric + H_contact` with a modified **Miyazawa–Jernigan** contact potential |
| qubits | dense encoding `2(N−1)`; a 7-mer is ~15–18 qubits; hardware threshold stated as **20 qubits**, depth 60–80 |
| largest run | **N = 10** residues, 135 sequences |
| structural accuracy | **none reported anywhere** — no Cα-RMSD, no backbone RMSD, no comparison to an experimental structure. Validation is energy-convergence plots and bitstring distributions |
| the "30% faster" claim | **convergence on the objective**, against a standard expectation-value VQE |
| classical baselines | **none** — no exhaustive enumeration, no simulated annealing, no greedy |
| simulability | not discussed |

Two consequences.

**First, it is not in conflict with anything here, because it measures a different thing.** "CVaR
converges 30% faster" is a statement about the *objective*. On this instrument the arm that
optimises the deployed objective hardest is **second-worst on structure** (Q3), and
`rho(optimisation gained, structure gained) = +0.003 [−0.068, +0.072]` across 36 matched-budget
cells. A convergence result carries no information about the structural question, in either
direction.

**Second, and this is the observation worth keeping.** At 15–20 qubits the full statevector is
32,768 to 1,048,576 amplitudes, and a tetrahedral lattice at N = 10 has `4^9 ≈ 262,144`
configurations — **the ground state is obtainable by brute-force enumeration in well under a
second**, so the VQE is competing against `min()` over a quarter of a million integers, and no
classical baseline is run. That is precisely the situation §3.1 measures for our own circuit, found
independently in a published package.

## 3.4 A fifth corroboration, and it is a field review rather than a group

`s21/BRIEF.md` asks for more external groups reproducing this programme's findings. The strongest
one is not in the quantum literature at all.

*Peptide–protein docking: from physics-based models to generative intelligence* (Chem. Commun.
62(28):7235, 2026; PMC13010362) states the position directly: **current peptide docking approaches
are often able to sample near-native binding modes, and correctly identifying those near-native
modes from decoys remains the difficulty**, because of the complexity of the binding energy
landscape. And a 2026 docking benchmark (arXiv:2605.03707) reports that on **8,597 systems** the
scoring function ranked the near-native pose top in **4,350 (50.6%)** and failed in **4,247
(49.4%)**.

That is `search-saturates-discrimination-binds` and `nothing-ranks-within-the-pool`, stated as
settled background by a **review** rather than by one group — which is a stronger form of
corroboration than a single reproduction, because it means the field has already priced it.

**Label it carefully.** Docking is a different task from de-novo peptide folding, the pools are
built differently, and their "near-native" band is not our in-band. What transfers is the **shape**
of the finding — sampling solved, selection not — and not any number. — **SUPPORTED** as an
external corroboration of the shape; **NOT** evidence about our magnitudes.

> **The published quantum peptide-folding literature reports energies and convergence; it does not
> report structure.** Of the four external efforts now examined, one (arXiv:2609.02113) reports
> Cα-RMSD — and the moment it does, it reports the energy-ranking failure this programme has
> measured six ways. The others report ground-state energies against Hamiltonians whose ground
> states are classically enumerable. **This programme's central negative is invisible to a
> literature that never measures the structure.** — recorded as an observation about the field,
> not a claim about any paper.

---

# 4. PRIORITY 4 — TWO RADICAL ALTERNATIVES

Judged by expected information per unit compute. Each carries a mechanism, a first experiment, a
matched control and a falsifier.

## 4.1 R1 — RUN, AND **DEAD**. A clean falsification of my own proposal, n = 126

`s21/d_divsel.py` → `d_divsel.json` + `.COMPLETE`. I did not leave R1 as a recommendation: it was
cheap, so I ran it. Basis **POINT CLOUD** throughout; the incumbent anchor is 3.0483, which
`d_cvarop` reproduced bit-exactly, so the comparison is against the production number and not
against a re-implementation. Fold-clustered CIs quoted.

**The pre-registered falsifier**, written into the module before the run: *dead unless `divmax`
beats **both** the incumbent top-75 **and** the random-in-band control by more than the MDE with a
CI excluding zero, **and** the sign control `divmin` loses.*

| band | `divmax` | vs top-75 | vs `rand_in_band` |
|---|---|---|---|
| top 30% | 3.0587 | +0.010 [−0.041, +0.057] | −0.015 [−0.044, +0.008] |
| top 50% | 3.0911 | +0.043 [−0.020, +0.124] | −0.046 [−0.115, +0.038] |
| all 500 | 3.7795 | **+0.731 [+0.571, +0.894]** | **+0.347 [+0.241, +0.445]** |

**All three bands DEAD. The falsifier fires exactly as written.** — **REFUTED**

**And the sign control goes the wrong way, which is the informative part.** At the unbanded pool,
`divmin` — *minimising* diversity — **beats** its matched random control by **−0.214 [−0.397,
−0.040], 66W/60L**. Minimising set diversity helps. So diversity is not the lever, my inferred
mechanism for §1.2d was wrong, and it is retracted there. The corrected mechanism is **error
coherence, not spread**: averaging cancels i.i.d. error and preserves systematic error, maximising
spread selects **outliers** rather than decorrelating errors (which is why `divmax` collapses to
3.78 on the unbanded pool), and minimising spread is outlier avoidance —
`consensus-is-outlier-avoidance` reproduced on a new operator.

**One thing survives the wreckage.** The score band is load-bearing: a random 75 drawn from the top
50% loses to the top-75 by **+0.0883 [+0.0113, +0.1631]**, CI excluding zero. The distogram's
ordering does real work *inside* the band — the fourth independent sighting, beside §1.2d's
`tail_member` −0.904, that the distogram orders in-band and nothing else in this project does.

## 4.2 R1 as originally proposed — kept unedited, because the pre-registration must not move

**Mechanism, and it comes out of §1.2d rather than out of the air.** The terminal operator is a
coordinate average over a selected set, and memory prices *where you average* at 1.024 Å against
*what you rank with* at 0.171 Å. §1.2d shows the averaging operator's value is destroyed by a
low-diversity set: Legacy, which selects compact pool-typical geometry, **helps** a single-member
readout by −0.41 and **hurts** the averaging readout by +0.33, against the same matched controls.
Every selector in this project ranks candidates independently and takes a top-k. **Nothing has ever
selected the set as a set.**

**First experiment (cheap, pool-restricted, no optimisation, ~15 minutes at n = 126).** Replace
"top-75 by distogram score" with "75 chosen to maximise a diversity criterion subject to a score
band" — a greedy determinantal / k-medoids-style pick over the pairwise-RMSD matrix `P` that the
instrument already computes once per target. Score band from the pool's own score distribution;
**no native information anywhere.**

**Matched controls, all three mandatory.** (i) top-75 by score — the incumbent, 3.0483 point cloud;
(ii) a **random** 75 from within the same score band, which is the operative control because a
diverse set and a random set are easy to confuse; (iii) 75 chosen to *minimise* diversity in the
same band — the sign control, which must go the other way.

**Falsifier.** If the diversity-selected set does not beat *both* the incumbent and the
random-in-band control by more than the MDE with a CI excluding zero at n = 126, R1 is dead.
Additionally dead if the minimise-diversity control does **not** lose, since that would mean the
criterion is not doing what it says.

**Prior against, stated first.** Memory `decorrelated-errors-exist-but-are-unusable` shows fusion
gains going as the *square* of the weaker channel's skill. That result is about fusing two
*predictors*; this is about choosing an averaging *set*, a different operator — but it is the
reason to expect a small effect, and R1 should be killed quickly if the first CI spans zero.

**Expected information per unit compute: high.** It reuses `P`, needs no energy evaluation, no
AMBER, no circuit, and it tests a lever nobody in this project has pulled.

## 4.3 R2 — PUT THE READOUT IN THE MATRIX. **OPEN, but its motivation is weakened by §4.1**

**Status first.** R2 as written below gates the readout on the tail's *spread*. §4.1 shows spread
is not the driver, so the gate as specified is the wrong variable and I am **not** recommending it
in that form. What survives is the premise — readout choice is a first-class variable with a 0.7 Å
swing on Legacy, larger than most of the levers this sprint is chasing — and the requirement that
any gate be an **error-coherence** proxy rather than a spread proxy. I do not have a native-free
coherence proxy, so R2 is recorded **OPEN and unmotivated** rather than promoted. The original text
is kept unedited below because a pre-registration must not move after data.


**Mechanism.** §1.2d makes readout choice a first-class variable with an effect size (**0.7 Å**
swing on Legacy) larger than most of the levers the sprint is chasing. The obvious construction is
per-target: use `tail_member`-like readouts where the tail is tight and `tail_avg` where it is
diverse. The gate must be **native-free**, and the natural one is the tail's own pairwise-RMSD
spread — already computed.

**First experiment.** On the same pool, for each target compute the tail's mean pairwise RMSD and
regress `(tail_avg − tail_member)` on it. Then a **fixed-global** readout vs a
**spread-gated** readout, decided by a threshold fitted on `dev_set(24)` and applied unchanged to
the 126.

**Matched control.** The fixed-global best readout (which is `tail_avg`), plus a **permuted-gate**
control that applies the same number of switches at random targets. The permuted gate is the arm
that says whether the *gate* has skill or the *switching* does.

**Falsifier.** If the spread-gated readout does not beat the best fixed-global readout by more than
the MDE with a CI excluding zero — or if the permuted-gate control matches it — R2 is dead. Also
dead if the threshold has to be re-picked on the 126 to work, which is a hyperparameter chosen on
the tuning instrument (§9).

**Prior against.** `score-axis-does-not-transfer` (11 ways of re-consuming the prior, none survived
dev) and `in-band-ordering-is-per-target` (per-target signals reach 0.986 within a target and 0.600
across) both say per-target gating is where this programme goes to lose. That is exactly why the
permuted-gate control is mandatory rather than optional.

## 4.35 R1's generalisation, WITHDRAWN before it could be quoted

I gave the coordinator a generalisation of T2 for the next sprint: *a candidate set generated from a
shared parametric prior has coherent errors and averages badly, while an independently retrieved set
has incoherent errors and averages well — a property of the GENERATOR, not the scorer or the set
size.* It was well-argued and it is **not supported by the one measurement that tested it** (§1.5f):
the latent's top-M averages against its matched random control by **−0.33 to −0.49 Å**, comparable
to or larger than the pool's **−0.377**, and it clears its own MDE at every M. **I withdraw it as a
next-sprint direction.** The caveat available to me — that the selectivities are not matched (pool
m = 75 of 500 = 15%, latent M = 64 of 8192 = 0.8%) — is one I did not register in advance, and I am
not rescuing a prediction with it.

## 4.4 NOT a third alternative — a steering correction that costs nothing

The brief allows at most two. R1 is dead and R2 is unmotivated by my own data, and I decline to
invent a third without a mechanism. But one measured fact from §1.2b is worth acting on and costs
nothing:

| α | ORACLE tail-average (the selector ceiling) | `disto` realised | **gap** |
|---|---|---|---|
| 0.05 | **1.644** | 3.075 | **1.431** |
| 0.15 | 1.963 | 3.048 (the incumbent) | 1.085 |
| 0.30 | 2.275 | 3.072 | 0.797 |

`disto`'s realised value is **flat in α** — 3.075 / 3.048 / 3.072, a spread of 0.027, well inside
any of these comparisons' own MDEs. The **ceiling is not flat**: it improves by 0.63 Å as α falls
from 0.30 to 0.05. So **the headroom available to a better in-band ranker is concentrated at small
α**, and it is nearly twice as large at α = 0.05 as at α = 0.30.

Concretely: **any ranker improvement should be scored at α = 0.05 as well as at the deployed
0.15**, because a ranker that is only slightly better will show almost nothing at α = 0.15 and up
to 1.43 Å of available room at α = 0.05. This is the ceiling-side counterpart of project memory's
`operator-consumes-set-mean` ("m\* shrinks 500 → 75 → 20 → 3–5 as the objective improves") — that
finding says the optimal set size shrinks as the ranker improves; this puts a number on what the
shrunken set is worth. Adding the α = 0.05 column to any selection table costs one line.

**And note what it is not.** 1.644 Å is an **ORACLE selector ceiling on a POINT-CLOUD basis**. It
is not achieved, not achievable without the native, and not comparable to the 3.204 built-chain
incumbent. It bounds a direction; it is not a target.

## 4.5 What I am NOT recommending

**Not a bigger or cleverer ansatz.** §3.1 closes it: the register admits no quantum resource at any
depth, and enlarging the register means leaving 9–16-residue peptides.

**Not more search.** Six instruments say search is saturated; §1.2's ceilings say the same from the
selection side (ORACLE tail-average at α = 0.05 is **1.644 Å** while `disto` achieves 3.075 —
the headroom is entirely in *which* candidates, not in *finding* more).

**Not another normalisation of Legacy + AMBER.** §1.2d says the readout swings the answer by 0.7 Å
on Legacy alone; until the readout is fixed and stated, a normalisation sweep is measuring the
wrong variable.

---

# 5. FURTHER AUDIT FINDINGS, produced after §2 was written

## 5.1 D13 — THE BAR ITSELF: "MDE = 0.084 Å" IS NOT A PROPERTY OF THE INSTRUMENT

`s21/d_mde.py` → `d_mde.json` + `.COMPLETE`. No structure rebuilt, no energy evaluated: the paired
differences already exist on disk.

`s21/BRIEF.md` §1 — and, identically, `s18/BRIEF.md` §6, `s19/BRIEF.md` and `s20/BRIEF.md` — states
*"MDE at 80% power = 0.084 Å. A null below that is uninformative, not negative."* Four sprints have
quoted it as a universal bar: *"5% of the instrument's 0.084 Å MDE"*, *"a null below the MDE"*,
*"a quarter of the MDE"*.

**Derive the operator before interpreting its statistic.** For a paired comparison at two-sided
α = 0.05 and 80% power,

    MDE = (z_0.975 + z_0.80) · SE  =  2.8016 · sd(paired differences) / sqrt(n)

`sd(paired differences)` belongs to the **comparison**, not to the instrument. The paired sd that
would make 0.084 Å correct at n = 126 is **0.3366 Å**.

**Measured across 26 real comparisons: per-comparison MDE / 0.084 spans 0.11× to 9.54× — a factor
of 84 — and only 12% fall within a factor 1.5 of the constant.** My pre-registered falsifier (*D13
is REFUTED if they cluster within ~1.5×*) does not fire.

**It errs in both directions, and both failure modes are in the table.**

| comparison | n | effect | sd(paired) | SE | **its own MDE** | ×0.084 |
|---|---|---|---|---|---|---|
| AMBER on the built chain (`rmsd_full − rmsd_arm`) | 126 | +0.0207 | 0.0385 | 0.0034 | **0.0096** | **0.11×** |
| projection λ=0.3 vs point cloud | 126 | +0.1664 | 0.2046 | 0.0182 | 0.0511 | 0.61× |
| `[R1] divmax(0.3) − rand_in_band` | 126 | −0.0151 | 0.2378 | 0.0212 | 0.0593 | 0.71× |
| **my own** `[D1] LEG pooled emb−θ` | 10 | +0.0448 | 0.2729 | 0.0863 | **0.2418** | 2.88× |
| `[D2] tail_avg(0.15) − argmin \| legacy` | 126 | −0.7354 | 2.1520 | 0.1917 | 0.5371 | 6.39× |
| **`[s20] encoding AMBc/spsa`** | 10 | −0.6492 | 0.9045 | 0.2860 | **0.8014** | **9.54×** |

* **Too large for a low-variance comparison.** The deployed AMBER repair is a **6.09 SE** effect at
  power 1.000. The brief's §3 calls it *"a quarter of the MDE"*, which invites reading one of the
  most precisely measured quantities in the programme as negligible. Sprint 20 labelled it
  ESTABLISHED and was right to; the *bar* is what misleads.
* **Too small for a high-variance comparison, and this one is live.** Sprint 20's `AMBc/spsa`
  encoding cell reported **−0.649 [−1.205, −0.165], "sig"** — while the panel's own 80%-power MDE
  is **0.8014 Å, larger than the effect it detected.** `|d|/SE = 2.27`, post-hoc power **0.62**,
  Type-M magnitude exaggeration conditional on significance **1.27×** (Type-S ≈ 0, so the sign is
  safe). The `LEG/spsa` cell is worse: power **0.18**, Type-M **2.39×**. That is a **third**
  independent reason the Sprint-20 encoding panel is fragile, beside the Nelder duplicate (§1.3a)
  and the step-count confound (§1.3c).
* **And it corrects me first.** My own LEG null has MDE **0.242 Å**, not 0.084. The **NOT MEASURED**
  label in D8 is right, but its bar was wrong: **my LEG null is uninformative below 0.242 Å**, and
  §1.3c should be read with that number. I have weakened my own conclusion rather than keep the
  flattering bar.

**The fix costs one line per comparison: report SE beside every mean and quote the MDE the
comparison's own paired sd implies.** The programme derived this correctly once —
`s16/review_FINDINGS.md` Q4.1 computes *"SE = 0.051 Å, MDE at 80% power ≈ 0.16 Å, power at a true
0.05 Å effect ≈ 9%"* for its own panel, and `s16/LEDGER.md:678` records it the same way. **The
constant is what generalised into `s18/BRIEF.md`; the method is what should have.**
— **ESTABLISHED**, methodological.

---

## 5.2 D14 — LEGACY'S IN-BAND RANK SKILL IS THE DISTOGRAM'S

`s21/d_partial.py` → `d_partial.json` + `.COMPLETE`. **Cost: zero.** Every coefficient was already
in Workstream A's live artefact; nothing was rebuilt.

Workstream A's `tailprice` report prints, as a headline distributional fact:

    ORACLE rank skill against TRUE RMSD:   disto +0.6332 | legacy +0.3348 | amber -0.0191

which reads as *genuine Legacy carries about half the distogram's in-band rank skill*. Three lines
above, the **same artefact** records `rho_dis_leg = +0.6452` median — Legacy is strongly rank
correlated with **the selector already in production**. §9: derive the operator before interpreting
its statistic. A marginal correlation with the truth, for a score already correlated with the
deployed score, is not that score's contribution. The quantity that answers *does Legacy add
anything* is the partial, computed **per target** from that target's own three Spearman
coefficients and aggregated **fold-clustered**:

| quantity | mean | CI95 fold-clustered | median | n negative |
|---|---|---|---|---|
| `rho(disto, true)` | +0.6179 | [+0.5233, +0.7072] | +0.7417 | 3/42 |
| `rho(legacy, true)` **marginal** | **+0.3281** | [+0.3011, +0.3552] | +0.3074 | 8/42 |
| `rho(disto, legacy)` | +0.4784 | [+0.4226, +0.5240] | +0.6243 | 5/42 |
| **`rho(legacy, true \| disto)` PARTIAL** | **−0.0076** | **[−0.0616, +0.0681]** | −0.0332 | **25/42** |
| `rho(amber, true \| disto)` PARTIAL — **control** | −0.0290 | [−0.0593, −0.0036] | −0.0322 | 25/42 |

**Legacy retains −2% of its marginal rank skill.** My pre-registered falsifier — *REFUTED if the
partial keeps more than half the marginal (> +0.16) with a CI excluding zero* — does not fire.
AMBER's marginal is already ~0 and its partial stays ~0, so the partialling is not manufacturing
structure.

**Independent cross-check that needs no partial at all:** Workstream A's own tables show `d+l`
never beating `d` at any α or any readout. This is the *mechanism* for that, and it is the same
class of error as Q11 (every circuit-side landscape metric collapsed once target difficulty was
partialled out) and as the shared referent floor.

**Consequence for §4 of the brief:** the mandatory matrix's Legacy row must be read as a partial,
not a marginal, or Legacy will look like a contributor when it is a correlate.
— **ESTABLISHED**, n = 42 (Workstream A's artefact was at 42/126 when this was computed; the
completion flag records the n and says to re-run when it reaches 126).

**A label defect in the same live report, two words.** `tailprice.py` annotates that block
*"(negative = lower energy is better)"*. The sign is inverted: `disto`, the score that works, is
**+0.633** with 2/35 targets negative, so **positive** rho means correct ordering. As printed, the
legend says the arm that works is the arm that fails.

## 5.3 CROSS-VALIDATION BANKED: §1.2d is now a two-lane result

Workstream A rebuilt `tailprice.py` with medoid and random-member readouts and re-primaried on
arm-minus-matched-random. Its numbers and mine were produced by different code, on different
target counts, from different intermediate representations:

| comparison, α = 0.15 | me (`d_cvarop`, n = 126) | Workstream A (`tailprice`, n = 42) |
|---|---|---|
| `tail_member \| legacy` | −0.407 [−0.557, −0.259] | −0.406 [−0.682, −0.122] |
| `tail_medoid \| legacy` | +0.302 [+0.204, +0.404] | +0.305 [+0.097, +0.530] |
| `tail_member \| disto` | −0.904 [−1.078, −0.729] | −0.874 [−1.166, −0.572] |

**Legacy beats a matched random subset on a single-member readout and loses to it on the medoid,
in both lanes.** And the pipeline identity of §1.2a is now confirmed in **three** places — A's
`tail0.15|disto`, my `tail_avg0.15|disto`, and production's `rmsd_avg` agree with
`max|diff| = 0.000000` on all 42 targets A has completed.

**One more thing from A's table, banked because it settles a §4 row.** AMBER is **worse than its
matched-count random control** at essentially every readout and every α — averaged tail
**+1.614 / +0.998 / +0.485** at α = 0.01/0.05/0.15 with **0/5 folds negative**, medoid up to
**+2.008** — and its ORACLE rank skill against true RMSD is **−0.042 [−0.116, +0.032]** with 23/42
targets negative. **AMBER contributes nothing to selection over this pool and actively hurts at
small α.** That is L4c ("composites never beat their own first stage") reproduced on a tail
selector.

---

# 6. WHAT I LEAVE OPEN, AND WHAT WOULD CLOSE IT

| # | open item | what closes it |
|---|---|---|
| ~~O1~~ | **CLOSED (§1.3d).** D1 on AMBc: −0.295 [−0.636, +0.022] at exactly equal step counts, 4 seeds, `.COMPLETE`. Falsifier did not fire. | — |
| O1b | **Whether a real encoding effect below ~0.5 Å exists on AMBc.** The panel's own MDE is ~0.53 Å, so §1.3d cannot exclude an effect of the size Sprint 20 reported, only show that it does not replicate. | The AMBc cell at n ≳ 60, which is ~6× the AMBER cost of the run I could afford tonight. Or, far cheaper, **fix the gauge bug (D16) first** — an embedded arm with a tangent-projected update is a different experiment and the current one does not bound it. |
| ~~O2~~ | **CLOSED (§1.3f), and not in my favour.** The step response is significantly signed in **both** directions (LEG +0.0912, DIST −0.0865, both CIs excluding zero) and predicts the encoding effect to **+0.004 Å** on DIST. The confound is fully explained; my stated mechanism (*"fewer steps helps here"*) is **REFUTED**. | — |
| O3 | **The unconverged CVaR law.** §1.2e's EXACT limit assumes a converged law. `tail_member` prices the uniform worst case; a real partially-trained law sits somewhere between `tail_member` and `argmin`, and nothing measures where. | An arm that reads out `argmin over the law's SAMPLES` at a realistic shot count, against `tail_member` and `argmin` as the two ends. Cheap. |
| O4 | **Whether AMBER's readout sign also flips.** §1.2d is measured on `disto` and `legacy`; genuine AMBER over 126 × 500 continuous single points was deliberately not spent on an operator question. Given L2c, AMBER (which prefers *expanded* geometry) may show the opposite pattern to Legacy — which would be a strong test of the compactness mechanism. | The same four readouts with `H = amber`, ~63,000 single points. Workstream A's corrected `tailprice.py` now carries `tailmem`/`tailmed` arms and may deliver it. |
| O5 | **The exactly-zero-SPSA-gradient rate**, flagged OPEN by my predecessor and still unmeasured: common random numbers plus a discrete sample space plus a tail selector can give `fp == fm` exactly. One line (`count fp == fm`). | Any lane running SPSA adding the counter. |
| ~~O6~~ | **CLOSED (§3.3).** QuPepFold read for methods: tetrahedral lattice + MJ contacts, 15–20 qubits, N ≤ 10, **no structural accuracy metric anywhere**, no classical baseline, "30% faster" is objective convergence. | — |
| ~~O7~~ | **CLOSED (§1.5a).** Built by the coordinator, extended by me to `n = 14–16`, **complete at n = 126**: +1.8149 [+1.6129, +2.0594], **0 of 122 targets reversed**. | — |
| O9 | **The cross-`n` magnitude of the discrimination gap.** Confounded by min-of-N on the ORACLE arm (a min over `2^16` versus `2^13`), declared in advance and labelled NOT MEASURED whichever way it landed. | `s21/latentrank.py`'s `oracle_8192` column — a min over a fixed seeded 8,192-configuration subset, identical in size at every `n`, coinciding with `latent_oracle` by construction at `n ≤ 13`. Running. |
| O10 | **Where in the objective's own ranking the ORACLE-best sits** — the question that decides whether a better *readout* can reach it at all. Registered prediction in `PREREG_D.md` D9: **0.1st to 5th percentile**, with a falsifier on both sides. **And a hazard raised before its data exists**: the top-M ceiling is a **min over M**, so it falls with M for reasons unrelated to ordering skill; the matched control is a **min-of-M over a random window of the same M**, and the vertical gap between the two curves is the ordering skill in the currency a reranker spends. | `s21/latentrank.py` completing, with the random-window control added. |
| O8 | **The gauge fix (D16).** `_EmbField` injects FD noise into a pure gauge direction and anneals its own step by ~2×. Every embedded arm this programme writes will inherit it. | One line in `s20/qb2_run.py`: project the update onto the tangent, or renormalise `u` per iteration. Then D1's AMBc cell is worth rerunning. |

---

# 7. THE STANDING FALSIFIERS, RESTATED SO THEY CANNOT MOVE

* **D1.** If, at matched realised step count, the embedding still wins by more than the MDE with a
  CI excluding zero **in both matching directions**, the encoding lever is not a step-count effect.
  **I withdraw D1 and the lever is SUPPORTED.** (Registered in `PREREG_D.md` before any grid ran.)
  **Outcome: did NOT fire in any cell.** Final numbers: `AMBc/spsa` −0.295 [−0.636, +0.022] at
  256/256 steps (n = 10); pooled over the **powered** n = 60 replication **−0.0061 [−0.0448,
  +0.0329], 30W/30L, sign p = 1.0000**, own MDE 0.055; `LEG/adam` −0.039 [−0.090, +0.010] and
  `DIST/adam` +0.004 [−0.032, +0.039] at matched steps. D1 stands; the lever is **NOT SUPPORTED**.
  **My own mechanism, however, is REFUTED in its stated form** (§1.3f): the step response is
  significantly signed in *both* directions, so *"fewer steps helps here"* is false as a general
  claim — what is true is that the encoding effect equals *minus the objective's own step
  response*, whatever its sign, and that prediction lands within 0.004 Å on DIST.
* **D2b.** If the zero-information `rand` energy's `tail(0.15) − argmin` gap were **less than half**
  the physics gaps, the endpoint would not be dominated by averaging and D2b would be REFUTED. It is
  **1.41×** the weakest physics gap. Not fired.
* **D3.** If the λ-continuation beats its **`toward_member`** displacement-matched control by more
  than the MDE with a CI excluding zero, the continuation carries direction information and D3 is
  REFUTED. Registered before Workstream C's artefact existed. **Beating `rand_iso` does not count**,
  because AMBER already does that by −0.378.
* **D8 (registered before the `n = 14–16` enumeration).** P1 — the gap is larger on `n ≥ 14`, by
  min-of-N on the ORACLE arm: **direction right (+0.252 / +0.301, inside my +0.2–0.6 band), CI on
  the difference spans zero → NOT MEASURED as pre-committed, and the MECHANISM is REFUTED**
  (`corr(n, latent_oracle) = +0.297`, the wrong way). P2 — qualitative conclusion unchanged:
  **CONFIRMED**, 0/51 new and 0/122 merged. P3 — exhaustive search still does not beat the pool on
  the targets where the sampler was genuinely searching: **CONFIRMED**, +0.1738 [−0.0079, +0.3715].
* **D10 (registered before `s21/d_lrank.py` produced a number).** T1 — medoid/avg beat the deployed
  argmin inside the top-M: **SUPPORTED.** T2 — the avg-vs-matched-random gain is *smaller* on the
  latent than on the pool and may not clear the MDE: **REFUTED on both halves** (−0.478 at M = 64,
  −0.327 at M = 512 against the pool's −0.377, clearing its own MDE at every M). T3 — no native-free
  readout reaches the ORACLE top-M ceiling: **SUPPORTED** (1.05 Å unclosed at M = 512).
* **D11 (registered before any AUC was read).** The bimodal containment route is DEMONSTRATED only
  if the best native-free predictor's folded |AUC| exceeds the **95th percentile of a permuted
  best-of-K maximum** (≈ 0.690 at this n and split) and its own CI excludes 0.5. **Anything in
  0.55–0.65 is NOT MEASURED and the route NOT DEMONSTRATED — never "promising".** Simulated at this
  experiment's own n: pure noise gives a best-of-10 |AUC| whose *median* is **0.625**.
* **D9 (registered before `s21/latentrank.py` produced a number).** The ORACLE-best's percentile in
  the objective's own ranking lands **between the 0.1st and 5th**; worse than the 20th or better
  than the 0.01st **refutes** it. And the top-M ceiling **must** be plotted against a min-of-M over
  a **random window of the same M**, or a curve that falls only because 512 draws beat 1 draw will
  be read as ordering skill.
* **D4.** No falsifier is claimed for a literature search. The reachability verdict in §3.1 is a
  measurement and its falsifier is arithmetic: produce a register size or a circuit family used on
  this instrument whose complete law cannot be tabulated.
