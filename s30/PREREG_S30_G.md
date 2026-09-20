# PREREG — S30 lane G

Two questions. Both falsifiers below are registered **before the first number exists**; this file is
committed before any measurement script is run. Contract rules 2, 3, 4, 5, 6, 7, 9, 12, 26 apply.

Author: lane G. Date stamped at commit (`date` run in the same command).

---

## Q1 — WHERE DOES THE k=30 AMBER-RELAX GAIN LIVE?

### The literature claim under test

Lane L found that the averaging-artifacts literature asserts *"averaging artifacts become more
pronounced when members of the ensemble are more divergent"* and **never measures it**. Our pool's
divergence is observable per target and the relax measurement already exists, so the claim is
testable here at zero compute.

### The artefact (verified before planning, contract rule 13)

```
ls  s16/results/repair_A.json        1,192,035 b   2026-09-06 07:33   126 per_target rows
ls  s16/results/repair_report.json      85,841 b   2026-09-06 09:10
git log --all --oneline -- <both>    a15406c8 checkpoint: pre-deep-dive research baseline
```

`repair_report.json.settings.k30_full.vs_proj_ungated` reads
`mean_diff -0.0221, se 0.0070, n 126, n_better 67, n_worse 59` — which is exactly the
`-0.022 [-0.036, -0.009]` the project's memory records for this arm. **The arm is identified and
nothing is recomputed.**

Per-target effect:  `d_t = rmsd(f0_k30_full)_t - rmsd(proj)_t`.  Negative = the relax helps.
Baseline `proj` is `I.project`'s `fit_ca`, i.e. the arm's own conservative baseline.

### Divergence variables — native-free, all measured on the set that is averaged (the shipped top-75)

| tag | definition | source |
|---|---|---|
| **DISP_rmsd** (primary) | mean pairwise CA-RMSD among the shipped top-75 | computed from `s12.instrument` |
| DISP_S (secondary) | RMS deviation of members from the set centroid | lane X `s30_X_typicalgood.json.pool[].S_spread` |
| DISP_rg (tertiary) | sd of Rg over the top-75 | lane F `s30_F_score.json.per_target[].F4_top75_rg_sd` |

DISP_rg is named third on purpose: lane F measured that pool Rg dispersion does **not** predict
filter failure (-0.128, CI includes zero), and that is what made lane L lower its own odds.

### Strata

`FAIL18` / other 108, from `s12/instrument.py:49`. **Declared limitation, stated before the
numbers:** FAIL18 is defined by the *filter's own in-band recall* (S30-L23), not by `d_t`, so this
is not the S30-L8 failure mode — but it *is* correlated with the production RMSD level, so the
FAIL18 row is reported as a **descriptive split, not a clean causal stratum**, and
fold 0 contains no FAIL18 target so its CI draws from four clusters and is wide by construction.
An ORACLE-defined tail (worst-18 by pool mean) is reported beside it and **labelled ORACLE**.

### POWER, computed and registered BEFORE the split (contract rule 3)

`sd(d_t) = se * sqrt(126) = 0.0070 * 11.225 = 0.0786 A`.

```
half-split difference    SE ~ 0.0786*sqrt(2/63) = 0.0140    MDE = 2.8016*SE = 0.0392 A
FAIL18 vs 108            SE ~ 0.0786*sqrt(1/18+1/108) = 0.0200    MDE = 0.0561 A
whole-sample effect                                             -0.0221 A
```

> **A half-split cannot reach its own MDE unless MORE than 100% of a 0.022 A effect sits in one
> half.** This arm is underpowered for localisation by construction. That is registered here, in
> advance, so that a null cannot be read as evidence of uniformity and a positive cannot be read
> without the multiplicity it costs.

Consequence for the reading: the **continuous** statistic is primary (it uses all 126 targets and
does not throw away the ordering), the split is descriptive, and concentration is scored against a
**uniform-effect null** (contract rule 5), never against a raw threshold.

### Falsifiers

- **F-G1 (primary, continuous).** Fold-clustered Spearman `rho(d_t, DISP_rmsd)` across 126
  targets, against a **fold-preserving permutation null** (10,000 draws: permute the dispersion
  labels within fold). **FIRES** if `rho <= -0.7 * (2.8016 * sigma_null)` — i.e. the relax gain
  grows with pool divergence — with a fold-clustered CI excluding zero.
- **F-G1b (split, descriptive).** `mean d_t` on the high-dispersion half minus the low half.
  **FIRES** only at `<= -0.7 * 0.0392 = -0.0274 A` with a fold-clustered CI excluding zero.
- **F-G1c (stratum).** `mean d_t` on FAIL18 minus the 108. **FIRES** at `<= -0.7 * 0.0561 =
  -0.0393 A` with a fold CI excluding zero.
- **Concentration** is reported as the observed high-minus-low gap's percentile under a
  **uniform-effect null** (same marginal distribution of `d_t`, dispersion labels permuted), plus
  median, mean/sd and win-rate **as one verdict** with the note that a nonlinear statistic has three
  different nulls (S30-L4). `d_t` is linear in RMSD so the mean is the right primary here; the
  median and W/L are reported beside it, not instead of it.
- **Mechanism sub-check (free, no extra bar).** Does divergence predict the *artifact* at all?
  `rho(DISP_rmsd, contraction of the raw average)` and `rho(DISP_rmsd, raw_avg_rama_ok)`, from the
  `raw_avg_geom` block already in `repair_A.json`. If divergence does not even predict the artifact
  on our instrument, the literature claim is untestable here rather than false, and that is the
  honest verdict.

### MY REGISTERED PRIOR, before looking — "IT DOES NOT CONCENTRATE", about 2 to 1 (i.e. ~1/3 that it does)

Lane L predicted "it concentrates" at 2 to 1, then lowered itself to roughly even after lane F's
Rg-dispersion null. I register **lower than lane L's revised odds**, for three reasons that are
mine and not lane L's:

1. The mechanism of this arm is **contraction repair**. Coordinate averaging contracts the backbone
   25.8% (`averaging-space-beats-the-objective`), and the contraction of a mean is set by the
   members' **common-mode alignment** and their number, not by their spread. Lane X measured
   `corr(S, B) = +0.0851` within target — spread is **orthogonal** to the bias the endpoint is made
   of, so there is no algebraic reason for a bias-repairing operator to scale with spread.
2. `pool-error-is-68-percent-common-mode`: the part of the pool's error that averaging leaves behind
   is the part that is *shared*, and shared error is by definition the low-divergence component.
3. The project's own record already shows this arm concentrates on **a few large winners**
   (n=50: top-10 share 1.163) — but "concentrated on the biggest winners" is a statement about the
   outcome and is **not** the same as "concentrated on the most divergent pools". Reading one as
   the other is exactly the S30-L8 error in a new place.

### Decision rule

If neither F-G1 nor F-G1b fires: **the literature's divergence claim is NOT SUPPORTED on our
instrument**, escape E2 is **uniform in pool divergence** and is therefore **not** a tail
intervention, and E2 closes as a route to the tail prize. Reported with the power note attached, so
"we could not localise it" is never dressed up as "it is uniform".

---

## Q2 — CAN A NATIVE-FREE CHANNEL WITH GENUINE GLOBAL REACH BE CONSTRUCTED AT ALL?

### The claim I am registering, before any measurement: THE CHIRALITY DICHOTOMY

**Setting.** A deployable single-structure channel is a map `S: R^{n x 3} -> R` invariant under
rotation and translation (a score must not depend on where the molecule sits).

> **THEOREM (G1).** For such an `S`, `S` is a function of the pairwise distance matrix `D` **iff**
> `S` is also invariant under reflection.
>
> *Proof.* (<=) Classical MDS: `G = -1/2 J D^2 J = X X^T` for the centred coordinates `X`, so `D`
> determines `X` up to right-multiplication by `Q in O(3)`. `O(3) = SO(3) x {+-I}`, so the residual
> ambiguity after rotation is **exactly reflection**; a reflection-invariant `S` is therefore
> constant on the fibre and descends to a function of `D`. (=>) `D` is manifestly
> reflection-invariant, so any function of it is. **QED**

**Corollary G1a — every candidate on the brief's list is a distance-map reading.** Long-range
contact topology, the radius-of-gyration *profile*, principal-axis / inertia-tensor structure
(asphericity, acylindricity, shape anisotropy), end-to-end and intermediate-separation distance
distributions, lane L's separation profile, contact order, excluded volume, packing density, and
**burial / SASA** are each **reflection-invariant** — a mirrored peptide has identical contacts,
identical Rg profile, identical inertia spectrum, identical separation profile and identical
solvent-accessible surface. By G1 every one of them is a function of `D`. They are not new
channels; they are **coordinate systems on the distance map**.

**Corollary G1b — the reference, not the functional, is where the information is.** A function of
`D` alone is a shape descriptor; to score *nativeness* it must be compared to an expectation. This
project has exactly three sources of an expectation, and the brief's three closed buckets are
exactly those three:

| reference | bucket | why it is closed |
|---|---|---|
| the predicted distogram (sequence channel) | (a) class M | bounded by theorem 2 |
| the pool | (b) consensus | measures typicality; anchor contrasts CONS -0.056, DMAP_CONS -0.053, POOLGO -0.050 |
| universal physics (AMBER, DOPE/LEG, Ramachandran, a hydrophobicity table) | — | **target-independent by construction**, so it cannot supply the per-target *sign*, which `in-band-ordering-is-per-target` says is the only leverage; and it is exactly what lane R's ladder matches rung-for-rung |

**Therefore the brief's three buckets are complete up to EXACTLY ONE family: CHIRAL functionals.**
A chiral functional distinguishes a structure from its mirror image, so by G1 it is **provably not a
function of `D`**, hence provably not a distogram re-reading; it needs no reference set, hence not a
pool-consensus term; and it is a non-separable integral over the whole chain, hence not a sum of
per-residue terms. **It is the only construction I can devise that is outside all three buckets by
proof rather than by assertion, and the project has never built one.**

Note this also *explains* lane R's split verdict rather than merely surviving it: every one of the
43 channels lane R tested is either a per-residue sum (blind by D1) or a `D`-functional with one of
the three references. There was no fourth kind in the library to test.

### What I will build — no new instrument (brief's constraint)

`s30/s30_G_chiral.py` **imports** `s30/s30_R_ladder.py` and reuses its `Sampler`, `make_ladder`,
`native_torsions`, `chain_path`, `rama_cnt`, `rg_of`, `ca_rmsd_to`, `spearman`,
`partial_spearman`, and its seed rule (`crc32("s30R|<pdb>|<seed>")`) **in the same call order**, so
the ladder is bit-identical to the one lane R scored. Aggregation uses `s30_R_agg.paired_ci` and
lane R's own `contrast_A_minus_ANCHOR` and `contrast_pref` definitions. Only the channel list is new.

Channels:

- **WRITHE** — the Gauss double integral (discrete Klenin-Langowski) over the CA chain. Chiral.
- **CHIRAL3** — mean signed volume `det[x_j-x_i, x_k-x_j, x_l-x_k]` over backbone quadruples at
  separation >= 1, normalised by `Rg^3`. Chiral; this is the coordinate-space (not torsion-space)
  handedness, so it is NOT the per-residue Ramachandran sign.
- **CHIRAL3_LONG** — the same restricted to `|i-l| >= 5`, i.e. the part that is *only* reachable by
  a lever arm and cannot be produced by local helicity.

**THE MATCHED CONTROL, IN THE OPERATOR'S OWN SPACE (contract rule 6).** For each chiral channel `X`
its **achiral twin** `|X|` is scored in the same run. `|writhe|` is reflection-invariant, hence by
G1 a function of `D`; `writhe` is not. Same functional, same scale, same ladder, differing in
**nothing but chirality**. The chiral content of the channel is `X` minus its twin, and if the twin
matches `X` the chirality is doing nothing. Lane R's **anchor control** and **pool-member control**
travel with every number as well.

**AUDIT, asserted in code, not claimed in prose.** Reflect every structure (`x -> x * diag(1,1,-1)`)
and assert: `D` unchanged to 1e-9, `|X|` unchanged to 1e-9, `X -> -X` to 1e-9. This proves the
construction is chiral rather than asserting it.

### Falsifier

- **F-G2.** A chiral channel's `contrast_A_minus_ANCHOR` (fold-clustered, lane R's own statistic)
  **FIRES** at `>= +0.10` — lane R's registered margin, the bar LEG_torsion missed at +0.024 — with
  a fold-clustered CI excluding zero **and** exceeding its own achiral twin by `>= 0.7 x` the
  paired MDE of that comparison.
- **Secondary (preference, the thing lane R showed does not survive):** `contrast_pref`
  (`pref_near - pref_pool`) with a fold CI excluding zero.
- The floor travels with every sentence: lane R's torsion rebuild is **0.347 A**, not 0.

### MY REGISTERED PREDICTION: F-G2 DOES NOT FIRE, about 4 to 1

Mechanism, registered so the null is interpretable: at `n = 9-16` a backbone's writhe is dominated
by **local helical handedness**, and lane R's rungs draw from the fold's leakage-safe Ramachandran
table, which carries that handedness — so the chiral signal is **matched by construction** in the
same way the local statistics are. The genuinely non-local part of writhe requires a chain long
enough to cross itself, and a 13-mer barely does. I expect `contrast_A_minus_ANCHOR` in
**0.00 to 0.05**, i.e. the neighbourhood LEG_torsion landed in, for a structurally similar reason.
`CHIRAL3_LONG` is the arm that could surprise me, because it is defined to exclude the local part.

### What the answer is if F-G2 does not fire

> Then the answer to the brief's question is **no**, and it is **a theorem rather than an
> enumeration**: by G1 every achiral native-free channel is a distance-map reading, the only three
> available references are the three closed buckets, and the single family that escapes the theorem
> is measured on this project's own instrument to be inside the noise at peptide length. That
> completes the sprint's ceiling argument instead of adding to it.

Stated as its own limitation: G1 is a statement about **single-structure** channels. It does not
bound set-referenced channels (which are bucket (b)) and it does not bound a channel with access to
a fourth reference this project does not have (e.g. an experimental observable). Both exclusions are
named here so the theorem is not over-quoted.
