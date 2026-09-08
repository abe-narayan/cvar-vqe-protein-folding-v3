# SPRINT 18 — PHYSICS PRE-REGISTRATION

**Written 2026-09-06, before any Sprint-18 physics number existed.** Binding. Nothing below is
adjusted after a result is seen; anything that changes is recorded as a change with its reason.

Workstream: PHYSICS (Legacy + AMBER). Phases 6, 8, 10 of the coordinator's brief.

---

## 0. THE HYPOTHESIS, STATED SO IT CAN FAIL

Sprint 17 (`s17/phys_FINDINGS.md` §3.2, ledger L29) measured on 126 targets × K = 500 identical
candidates:

* `leg_contact` (Miyazawa–Jernigan) carries **+0.080 [+0.032, +0.128]** of in-band partial
  Spearman with Cα-RMSD after the distance score is partialled out, against a **+0.010**
  matched-random null — real information;
* the same term **selects +0.1357 [+0.0433, +0.2327] Å WORSE than random in-band**, 52W/71L —
  it cannot be an argmin;
* its **global** Spearman with RMSD is **−0.176** — globally it is *anti*-correlated, i.e. the
  candidates it likes are the worse ones.

> **H-C (the workstream's headline hypothesis).** `leg_contact` fails as a *selector* but may
> succeed as a *term in the objective*. Adding it to the distance objective supplies the
> inter-residue structure that objective is missing (93% separable by Walsh spectrum against a
> truth at mean Pauli weight 3.72) and therefore moves the objective's optimum toward the native.

**I record my honest prior before running: I expect H-C to FAIL.** The global sign of
`leg_contact` is the wrong one (ρ = −0.176), the in-band information is a *bulk ordering* signal
and L29 already records that "a weak monotone ordering signal does not become an argmin", and an
objective's optimum is an argmin. The reason to run it anyway is that *term-in-objective* is a
genuinely different functional role that has never been tested, and Sprint 17's decisive negative
was measured on the ranking role only.

---

## 1. THE OBJECTIVE, AND WHERE ITS PIECES COME FROM

    E(theta; lam_c)  =  E_le1_dist(theta)  +  lam_c * s * E_contact(theta)

**`E_le1_dist` is NOT built here.** It is consumed through the *single shared* degree-1 ANOVA
definition owned by MATH (`s18/math_*.py`) via the adapter EXPERIMENT declared in
`s18/COORD_exp_to_math.md` (`s18/exp_anova.build(...) -> obj`, `obj.E_le1`, `obj.E_full`). If
MATH's definition supersedes the provisional one, **I re-run**; I never ship a second definition.
The full objective `obj.E_full` is the deployed distance objective, unchanged.

**`E_contact` IS built here and is genuine Legacy.** It is `core.energy.contact_term` evaluated on
the ideal-geometry backbone of (φ, ψ) — the *same* object `s16.energy_lib.legacy_components_of_windows`
computes, verified bit-for-bit against it before use (gate G0 below). `DEFAULT_WEIGHTS["contact"]
= 1.0` and **is never fitted**.

**Naming discipline.** EXPERIMENT's `lam` is a convex mix between `E_le1` and `E_full` and is a
*different quantity*. Mine is written `lam_c` everywhere.

### 1.1 NORMALISATION — declared exactly, and native-free

The two terms have unrelated units, so `lam_c` is meaningless without a stated scale match.

> **The scale is the per-target standard deviation of each term over the target's own
> ranker-neutral K = 500 BLOSUM62 retrieval pool** (`s12.instrument.load_univ`, retrieval order,
> the identical candidate set Sprint 17 used). Both terms are evaluated on the same 500
> ideal-geometry rebuilds. Then
>
>     s = sd_pool(E_le1_dist) / sd_pool(E_contact)
>
> so that `lam_c = 1` means **one pool-standard-deviation of contact energy trades against one
> pool-standard-deviation of distance objective**.

* This uses **no native information of any kind** — the pool is fixed by retrieval, and only
  candidate energies enter.
* The retrieval pool is used rather than the shipped top-75 because top-75 is selected *by the
  distance score*, which deflates that term's variance and would silently inflate `lam_c`. A
  top-75 sensitivity arm is reported beside the primary.
* Degenerate scale (`sd_pool(E_contact) < 1e-9`) excludes the target, counted and named.

### 1.2 THE SIGN — pre-registered before the run

`core.energy.contact_term` sums MJ-**corrected** pair energies weighted by a cosine switch, and
corrected MJ energies are **negative for favourable (hydrophobic) contacts**. Legacy is
*minimised* at weight +1.0. Therefore:

> **The physically motivated sign is `lam_c > 0`** — minimising rewards forming favourable
> contacts. It is also the sign of the +0.080 in-band partial correlation.
>
> **It is NOT the sign the global correlation favours.** ρ_global(`leg_contact`, RMSD) = −0.176
> says that over the whole pool, *lower* contact energy goes with *higher* RMSD; a global-fit
> reading would want `lam_c < 0`. **This contradiction is registered here, before the run, as the
> chief reason H-C may fail, and negative `lam_c` is carried only as a CONTROL — a negative
> result at `lam_c < 0` is not a rescue of the hypothesis and will not be reported as one.**

### 1.3 THE LADDER — fixed, and not extended

    lam_c  in  { -1, -0.5, 0, +0.5, +1, +2 }

Six arms, no more, whatever the result. `lam_c = 0` is degree-1 distance only. **`lam_c = +1` is
the single PRIMARY arm** — the physically motivated sign at unit variance match. Everything else
is shape. This is declared now so that a favourable non-primary arm cannot be promoted afterwards.

### 1.4 HOW `lam_c` IS (NOT) CHOSEN

**There is no native-free rule in this programme that selects `lam_c`, and I do not invent one.**
The whole ladder is reported. Any `lam_c` picked by looking at RMSD is labelled **ORACLE** and is
not a predictive result. The primary test is the single pre-declared `lam_c = +1`.

---

## 2. PHASE 6 — THE CRITICAL COMPARISON

On **identical starting candidates** (the shipped top-75 windows and their coordinate average,
`s14.avgspace.top75_windows`, and the ranker-neutral K = 500 pool for pool statistics), three
objectives:

    A   E_full        the deployed distance objective                (s17 refine_full's objective)
    B   E_le1         degree-1 distance only                          (= lam_c 0)
    C   E_le1 + lam_c * s * E_contact   for each lam_c in the ladder

Brief §5 forbids reporting these as one number. **Six axes, separately:**

| axis | statistic |
|---|---|
| **objective quality** | does optimisation lower the objective it is aimed at? (mean drop) |
| **objective alignment** | Spearman(objective, RMSD) over the K = 500 pool, GLOBAL and IN-BAND, and in-band selection vs matched random |
| **local refinement** | RMSD of the refined structure vs the coordinate-average start, paired |
| **near-native coverage** | P(pool has a sub-2 Å member), recall of sub-2 Å members |
| **candidate generation** | best member of the refined ensemble; ensemble RMSD after coordinate averaging |
| **ensemble diversity** | mean pairwise Cα-RMSD within the refined ensemble |
| **final RMSD** | full-chain Cα-RMSD of the emitted structure — **the primary endpoint** |

**PRIMARY ENDPOINT.** Final full-chain Cα-RMSD of the coordinate average of the ensemble refined
under `lam_c = +1`, against the same quantity at `lam_c = 0`, n = 126, target as the unit, paired
fold-clustered bootstrap CI, median and W/L reported beside the mean.

**SUCCESS (H-C SUPPORTED):** `lam_c = +1` beats `lam_c = 0` on the primary endpoint with a
fold-clustered CI excluding zero **and** beats both mandatory controls.

**FALSIFIER (H-C REFUTED):** `lam_c = +1` does not beat `lam_c = 0` with a fold-clustered CI
excluding zero. Then the branch closes: no ladder extension, no re-normalisation, no sign flip
sold as a result. Given the instrument's 0.084 Å minimum detectable effect at 80% power, a null
below that magnitude is reported as **uninformative**, not as negative.

**CONTROLS on every arm (both mandatory):**
1. **Zero-information**: the identical refinement with the contact term's MJ matrix replaced by a
   `stable_rng` **residue-label permutation** — same magnitude, same functional form, no sequence
   information. (This is the correct zero-information control for a *contact* term: shuffling the
   distogram controls the distance term, not this one.)
2. **Matched-random**: a random torsion move of the same realised magnitude as the arm's move.

---

## 3. PHASE 8 — THE CAUSAL DOWNSTREAM COMPARISON

Every candidate pool passes the **same** downstream evaluation on **identical structures**:

    no filter · leg_torsion · leg_contact · combined Legacy · AMBER single point ·
    AMBER restrained repair at k = 30 · Legacy -> AMBER · matched-random gate -> AMBER

Reported: mean RMSD, median, W/L, fold-clustered CI, near-native recall, best-member RMSD,
ensemble RMSD, diversity, Ramachandran (three-way), cis fraction and ω deviation, clash counts and
minimum heavy separation, Cα displacement, convergence count.

**Two standing traps, named in the brief and enforced here:**
* AMBER's validity statistics **do not** substitute for RMSD evidence. Every validity table is
  printed beside the RMSD table for the same structures and neither is quoted alone.
* Legacy's energy correlation **does not** substitute for selection evidence. Every correlation is
  printed beside its selection delta against a matched-random gate of the same count.

**AMBER discipline.** Convergence gate declared here, before use: *final energy finite and
≤ 1000 kcal/mol* (`core.amber.CONVERGE_MAX_KCAL`). Exclusion counts and the excluded PDB IDs are
reported with every gated arm. **A gated arm is compared to its own gated input**, never to the
ungated instrument-wide input. The rotated-frame null is reported with its **maximum**.

**Filter fractions** are pre-registered at f = 0.50 (the fraction Sprint 17 measured `leg_torsion`
at), with f = 0.25 as shape. Each gate is priced against a matched-random gate **of the same
count**, 3 `stable_rng` draws.

---

## 4. PHASE 10 — SPACING RESTORATION (secondary; may not crowd out 6 and 8)

The coordinate average contracts Cα–Cα spacing to **2.949 Å against an ideal 3.80 Å** and
Spearman(input spacing, cis fraction at `cafix`) = **−0.949**; the ideal-geometry projection costs
**+0.164 Å** to fix it. Arms, all native-free:

    S0  do nothing                                   (the 3.048 A average)
    S1  global uniform Ca rescale toward 3.80 A      (one parameter, no native)
    S2  segment-wise / local Ca spacing correction
    S3  minimal-displacement Ca correction (projection onto the spacing constraint set)
    S4  constrained reconstruction preserving the AVERAGED TORSIONS
    P   the incumbent ideal-geometry projection      (the +0.164 A reference)

**SUCCESS:** some arm reaches the projection's cis/convergence repair at **less than +0.164 Å** of
Cα cost, CI excluding zero against the projection. **Otherwise the tax is confirmed irreducible by
these means and the branch is reported closed.** No arm is promoted on validity alone.

---

## 5. GATES THAT MUST PASS BEFORE ANY SCIENTIFIC NUMBER IS QUOTED

* **G0 — the contact term reproduces genuine Legacy.** My torsion-space `E_contact(φ, ψ)` matches
  `s16.energy_lib.legacy_components_of_windows(seq, PHI, PSI)["contact"]` to < 1e-9 absolute on a
  pre-declared sample. If it does not, nothing downstream is reported.
* **G1 — the gradient is correct.** The finite-difference torsion gradient reproduces
  `s15.align_lib.fit`'s exact analytic gradient on the pure-distance objective to < 1e-5 relative.
* **G2 — the start reproduces the record.** The coordinate-average start reproduces
  3.0498 Å (n = 126) and its projection 3.213 Å, matching L30.
* **G3 — one definition of the degree-1 object.** `E_le1` is imported, never re-implemented, and
  the module records which implementation (MATH's or the provisional adapter) produced each
  artefact, by config hash.

## 6. ARTEFACT DISCIPLINE

Unique immutable directory per experiment under `s18/results/`; `stable_rng`/`stable_seed`
throughout; **explicit `complete` flag and row count on every write** (the hazard that let a 9-row
partial overwrite a 126-row artefact in Sprint 17); per-target checkpoint; persisted config hash
and the full `lam_c` ladder in every artefact. **The sealed 60-target benchmark is not read,
probed, or derived anywhere in this workstream.**
