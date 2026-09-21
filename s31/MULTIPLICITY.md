# SPRINT 31 — SPRINT-WIDE MULTIPLICITY REGISTER

**Owner: lane D.** Carried to S31 as a build item because in S30 each lane counted its own
comparisons and the sprint-wide total was only ever *a lower bound in the hundreds* — which is
not a number you can correct a p-value with.

## The rule

**Write to this file as comparisons are emitted, not at the end of the sprint.** A lane that
emits a family of comparisons appends one row per family *when the family is run*, whether or
not anything in it was positive. Comparisons that were run and thrown away count. Comparisons
inside a grid count as the size of the grid, not as one.

A row is:

| # | date time | lane | family | k (comparisons) | basis | what was claimed from it |

`k` is the number of *statistical comparisons the family could have produced a headline from*.
A 32-cell sweep scored on one endpoint is k = 32. A sweep scored on two bases is k = 64. A
single pre-registered contrast is k = 1.

## How the total is used

1. **Nothing in this register changes a pre-registered primary comparison.** A contrast named in
   a PREREG before any number was seen is k = 1 for its own purposes, and the register records
   it so the sprint total is honest.
2. **Exploratory families get the Bonferroni-adjusted MDE quoted beside the raw one.** With the
   project's fixed rule MDE = 2.8016 × SE (two-sided α = 0.05, power 0.80), a family of `k`
   comparisons needs roughly `z(1 − 0.05/2k) + z(0.80)` in place of `2.8016`; the adjusted
   multiplier is tabulated below so no lane has to derive it.
3. **A "best of k" result is an order statistic.** The standing S29 finding is that a per-target
   or per-cell minimum over k variants is *mostly best-of-k*; such a family must be priced
   against `ST.best_of_k_null` / `ST.best_of_k_within`, and the register row must say so.

### Adjusted MDE multiplier (two-sided, power 0.80)

Computed, not quoted: `norm.ppf(1 − 0.05/(2k)) + norm.ppf(0.80)`. The k = 1 row reproduces the
project's 2.8016 exactly (2.80159), which is the check that the formula is the right one.

| k | multiplier | vs 2.8016 |
|---|---|---|
| 1 | 2.802 | 1.00× |
| 2 | 3.083 | 1.10× |
| 4 | 3.339 | 1.19× |
| 8 | 3.576 | 1.28× |
| 16 | 3.797 | 1.36× |
| 32 | 4.004 | 1.43× |
| 64 | 4.201 | 1.50× |
| 128 | 4.388 | 1.57× |
| 256 | 4.567 | 1.63× |
| 512 | 4.738 | 1.69× |
| 912 | 4.876 | 1.74× |

(`k` = family size. Read: a 32-cell exploratory sweep needs **1.43× the single-comparison MDE**
before one of its cells is a result. S30's 32-cost meter sweep emitted 912 comparisons; at that
family size the multiplier is 1.74×. Bonferroni is conservative when the comparisons are
correlated — and in this project they are heavily correlated, because they are usually the same
126 targets scored different ways — so a lane that thinks the adjustment is too harsh should say
so **with the correlation measured**, not by preference.)

## Standing standards these comparisons are held to

- **MDE = 2.8016 × SE, per comparison.** MDE is *per comparison*, not a per-instrument constant;
  report SE beside every mean. The quoted 0.084 Å constant is wrong by up to 84× in both
  directions.
- **Below 0.7× MDE is not a result. 0.7–1.0× is NOT MEASURED.** Not "a trend", not "suggestive".
- **Fold-clustered CIs on the pinned folds.** Never regenerate folds or clusters.
- **ORACLE arms are labelled ORACLE / NOT DEPLOYABLE in the same sentence as their number.**
- **`s24.stats_lib.compare` is lower-is-better** (`d = a − b`, negative = a better). For a
  *preference rate* or any higher-is-better statistic its `verdict` string is **inverted**.
  **Quote the gate, never `.verdict`.** (Carried from S30.)
- **A control must match the operator's space**, and must be *plausible*, not uniform.
- **Any built-chain claim smaller than 0.0107 Å is inside the instrument's own reprojection
  noise** — see S31-L4 (defect D-B) and `s31/results/s31_D_projection_pin.json`.

---

## Register

| # | date time | lane | family | k | basis | claimed |
|---|-----------|------|--------|---|-------|---------|
| 1 | 2026-09-21 00:06 | D | D-B projection determinism: per-target pass-1 vs pass-2 bit-identity | 126 | built chain, 126 targets | determinism assertion only; no effect claimed |
| 2 | 2026-09-21 00:06 | D | D-B 1e-14 input perturbation vs unperturbed reprojection | 126 | built chain, 126 targets | instrument property (noise floor); no effect claimed |
| 3 | 2026-09-21 00:06 | D | D-B reprojection vs s29 canonical rows, and vs s27 DIS rows | 2 | built chain, 126 targets | agreement check; no effect claimed |

*(Rows 1–3 are instrument characterisation, not hypothesis tests: they carry no α budget because
nothing was accepted or rejected on them. They are registered anyway, because the rule is that
comparisons are registered when emitted and the decision about whether they spend α is made
here, in the open, rather than by whoever would benefit from the answer.)*
| 4 | 2026-09-20 23:57 | E | E1 direction: cos/corr(mu_hat, mu) levels + 2 matched controls, 2 pair masks | 28 | within-target cosine, 126 (106 long-range) targets | **1 pre-registered primary** (uncentered cos, all pairs, vs the permutation control) named in `PREREG_S31_E.md` §2 before any number; the other 27 are its robustness surface and carry no separate claim |
| 5 | 2026-09-20 23:57 | E | E2/E3 applied correction arms: 14 arms x {built chain, CA cloud} vs PROD | 28 | both bases, 126 targets | **5 pre-registered** (E2/E3 native-free; E2o/E3o/class-ceiling under AMENDMENT 1) + 3 mandatory matched controls; the ORACLE decomposition arms are bounds, not hypothesis tests |
| 6 | 2026-09-20 23:57 | E | E3 orthogonal SSE decomposition of y against mu and mu_hat | 18 | ORACLE energy/SSE shares, 126 targets | descriptive mechanism decomposition; **no effect accepted or rejected on it** -- it is the explanation of an effect S30 already measured |
| 7 | 2026-09-21 00:10 | D | S31-L8 branch experiment: B4 vs PROD, ORACLE_B4 vs PROD | 2 | built chain, 126 targets (the CA cloud is the INPUT and is identical across arms, so there is no second basis) | **both pre-registered in `s31/PREREG_S31_D_branch.md` before any arm was computed**, with the decision rule and the expected outcome (H0) fixed in advance. ORACLE_B4 −0.0938 Å, 2.92× MDE, **ORACLE / NOT DEPLOYABLE**. B4 −0.0055 Å, **0.44× MDE, NOT A RESULT**. Secondaries (conditioning, tail, FAIL18/108) are descriptive and named in the prereg. |
| 8 | 2026-09-21 00:10 | D | S31-L8 secondary strata: FAIL18 and the other 108 | 2 | built chain | pre-specified in the prereg §5.3; both reported (+0.0049 at 0.08×, −0.0072 at 0.69×), neither a result, **no stratum rescue attempted** |
| 9 | 2026-09-21 00:16 | B | S31-L10 B1 achirality: ff14SB phase table (F2), U(x) vs U(Rx) with a matched rotation control (F1), Legacy 11-term parity (F3) | 0 | potential-function property checks | **spends no alpha** -- these are properties of the force field, not comparisons against a null. Recorded because F1's registered 1e-6 threshold FIRED and the matched control inverted the reading; the row exists so that is not invisible. |
| 10 | 2026-09-21 00:16 | B | S31-L10 sweep parity: 19 chain/CA costs x {total, even, odd} preference contrast | 57 | built chain rungs, 126 targets | EXPLORATORY diagnostic of the S30 meter. No cost was accepted or rejected on it; it decomposes a published statistic. The reproduction of the published values to four decimals is a verification, not a claim. |
| 11 | 2026-09-21 00:16 | B | S31-L10 in-pool: 11 channels x 2 bands x {rho, rho|DIS, rho|RG} plus the G5 control contrasts | 78 | per-candidate CA point cloud, 126 targets | **G2, G3 and G5 are pre-registered primaries** in `s31/PREREG_S31_B.md` with falsifiers fixed before any number (k = 1 each for their own purposes). The remaining cells are exploratory and are reported at the adjusted multiplier. G2 fired, G3 did not, **G5 fired at 1.81-3.85x MDE**. |
| 12 | 2026-09-21 00:16 | B | S31-L10 candidate graph: C1 collapse (2) + C2 node statistics rho and rho\|medoid (9) | 11 | top-75 band, CA cloud, 126 targets | pre-check requested by the coordinator and run BEFORE anything was built on it. C1 did not fire, C2 fired. Nothing was deployed from it. |
| 13 | 2026-09-21 00:18 | F | **[renumbered in place from 9 at 00:20 -- collided with lane B's row 9, which is stamped 00:16; contract rule 18 renumbers the LATER-stamped entry, which is mine]**  S31-L11 F3 prefix-length order-statistic audit: F3a (chain+cloud), F3b zero-skill random m, F3c ORACLE-global + LFO m, F3d split-half transfer, F3e matched random-subset family (4 draws), F3f four native-free LFO m-rules, F3g effective-K | 12 | CA point cloud for the 126x128 matrix; built chain for F3a (basis stated on every row) | **pre-registered in `s31/PREREG_S31_F.md` s11 (commit de852bef) before any aggregate of the family existed**, with two bars written to fire. BOTH FIRED: the matched random family reaches 146% of the prefix gain (bar: >=80%), and the best fold-held-out native-free m-rule is -0.51x its own MDE (bar: >-0.7x). No positive claimed; the family DEFLATES a two-sprint-old number rather than asserting one. |
| 14 | 2026-09-21 00:21 | P | S31 lane P closed-form substitution: 6 arms (A quantum-OFF, B/C selection readout with VQE `p` / closed-form `p*`, D/E convex readout with `p` / `p*`, F shipped uniform ablation), 9 paired contrasts x {built chain, CA cloud} | 18 | **built chain is the primary basis** (canonical prod 3.2105 A); CA cloud carried as a second basis for comparability with lane A's cloud measurement and labelled as a different object on every row; 126 targets | **2 pre-registered primaries on the BUILT CHAIN ONLY** (`P1 = E-D`, `P2 = C-B`) named in `s31/PREREG_S31_P.md` §5 before any arm was computed, with the three-way gate (NULL / NOT MEASURED / MEASURED) and a 0.0107 A chain-floor override fixed in advance; k = 1 each for their own purposes. `S3 = D-A` is the built-chain cost of the SHIPPED quantum stage, promoted to a reported row by AMENDMENT 1 after the coordinator withdrew the +0.2260 A deficit (wrong readout). The other 6 chain contrasts and all 9 cloud contrasts are context and carry the adjusted multiplier. |
| 19 | 2026-09-21 00:58 | F | S31-L21 terminal operator: 3 operator arms x {chain, cloud}, 10 gates, the DISP mechanism contrast, the separation-band long/short contrasts, the coh contrast, and 4 tail strata x 21 statistics (descriptive), plus the UNREGISTERED ORACLE best-of-4 readout ceiling (declared as unregistered in S31-L21 s9b) | 29 | built chain primary, CA point cloud stated on every cloud row; the chain comparator is production RE-PROJECTED IN THE SAME JOB (3.2126), not S29's 3.2105 | **all falsifiers pre-registered in `s31/PREREG_S31_F.md` before the numbers existed (8a14edea, 49ee7c92, bc81f029, f5f69ba1)**. FOUR FIRED AGAINST THE LANE: MED +0.0688, AVG_RG +0.0475, AVG_SEP +0.4609, every gate worse in its registered direction, and the registered NEGATIVE dispersion mechanism measured POSITIVE. **No positive is claimed**; the family REFUTES this lane's own hypothesis. The tail table (4 strata x 21 statistics) is DESCRIPTIVE -- nothing was accepted or rejected on it -- so it spends no alpha and that judgement is recorded here in the open. |
| 21 | 2026-09-21 01:26 | F | S31-L11 s5d, the F3 arms on the BUILT CHAIN in ONE JOB: PREFIX(bestm128), ORACLE-global m, leave-fold-out m, and 2 matched-random-family draws, each vs M75, plus the 2 random draws vs PREFIX directly | 7 | built chain (the cloud is the selection basis for every arm and is reported beside it, with the per-arm cloud-to-chain price) | pre-registered in `s31/PREREG_S31_F.md` s11; **no positive claimed** -- the family carries S29-L30's transfer arms from the cloud onto the endpoint (-0.0044 at 0.19x, +0.0075 WORSE at 0.22x, both NOT A RESULT) and confirms the order-statistic bar on the chain at 140.6%. Draw 0 clears 1.0x MDE and draw 1 is at 0.93x; **both reported, claim rests on the draw MEAN**, contract rule 10. |

### Running total

**k = 759 comparisons registered** as of 2026-09-21 01:26 -- **DERIVED BY SUMMING ALL 24 ROWS**
(lane F's row 21 adds 7). Earlier totals below are left standing per rule 13.

**k = 706 comparisons registered** as of 2026-09-21 01:08 -- **DERIVED BY SUMMING ALL 20 ROWS**.
*(lane F's row 19 rose from 24 to 29 when it added the unregistered ORACLE best-of-4 arm; the
700 figure in this line's earlier version counted 24 and is corrected here.)*
*Correction, annotated in place (rule 13): the 611 line immediately below was mine and was WRONG
by 89. My regex required the comparison count to be the 4th pipe-delimited field, and lane B's
rows 11 and 12 contain literal `|` inside their descriptions (`rho|DIS`, `rho\|medoid`), so those
two rows -- 78 and 11 comparisons -- were silently dropped. A register whose own total does not
add up is worse than none, and that applies to my correction of it too.*

**k = 611 comparisons registered** as of 2026-09-21 01:00 -- **DERIVED BY SUMMING EVERY ROW'S
OWN COUNT** (17 rows), not by adding a delta to a previous line. The four lines below are
earlier totals left standing per rule 13; two of them disagree with each other because three
lanes edited this file inside four minutes. Lane F's row 19 adds 24.

**k = 478 comparisons registered** as of 2026-09-21 00:16
(the 332 below plus lane B's **146** in rows 9-12; the earlier 332 line is left standing
immediately after this one rather than overwritten)

**k = 490 comparisons registered** as of 2026-09-21 00:20
(lane B's 478 plus lane F's **12** in row 13. *Correction, annotated in place:* at 00:18
lane F overwrote the 332 line below while lane B was mid-edit and briefly published a
total of 344, which was 134 too low; lane B's 478 line above was never removed and the
332 line is restored here with its original wording.)

**k = 508 comparisons registered** as of 2026-09-21 00:21
(lane F's 490 plus lane P's **18** in row 14. The 490 line above is left standing
rather than overwritten -- `ledger-numbers-collide-under-parallel-lanes`: a running
total is a read-modify-write and this file has already had one collision tonight.)
| 15 | 2026-09-21 00:25 | C | S31-L5 C3 widening gate: 5 strata x 1 contrast x 2 bases (cloud, chain) + 1 random-18 null (20,000 draws, reused) + 2 diagnostics + 1 level regression + 2 fold-clustered paired contrasts | 13 | CA point cloud primary, built chain from s29 lane O rows (nothing reprojected); 126 targets | **F-C3a/F-C3b pre-registered in `s31/PREREG_S31_C.md` (commit 8ce5e1a0) before the first lane-C number**, with the -1.00 A bar and the null test fixed in advance; k = 1 for its own purposes. F-C3a FIRED at -1.1811 A. The 3 non-headline strata are pre-specified diagnostics, not rescue attempts, and `defn18` was registered as a definition-matched control BEFORE it turned out to equal FAIL18 exactly. |
| 16 | 2026-09-21 00:25 | C | S31-L14 set-matched ladder: 11 named arms + 21-point m-curve on 2 supports (42) + 8 fold-clustered contrasts vs the ORACLE argmin + 33 vs production + 21-value ridge path + 8 native-free/3 ORACLE cross terms | 76 | CA point cloud, 126 targets; the BUILT CHAIN block re-uses s29 lane O rows and reprojects nothing | **F-C1a, F-C1b, F-C1c, F-C1d pre-registered** (8ce5e1a0 / amendment 2e12e02c), k = 1 each. F-C1a REFUTED (-0.0101 A, 0.56x MDE), F-C1b CONFIRMED, F-C1c REFUTED (out-of-fold R2 0.021 vs a 5% bar, applied +0.0309 WORSE), F-C1d REFUTED (the affine ceiling is a dimension-counting tautology). The m-curve, the ridge path and the cross-term table are DESCRIPTIVE decompositions -- nothing is accepted or rejected on them -- and no per-target maximum is reported anywhere in the entry. |
| 17 | 2026-09-21 00:25 | C | S31-L15 candidate index: 6 maps x (1 rho + 1 ORACLE ceiling + 6 j-levels x 4 cell statistics + 7 bit MIs) = 222 descriptive, of which the 72-cell (map, j, rule) grid is alpha-spending; plus 5 fold-clustered contrasts for F-C2a | 77 | CA point cloud, 126 targets | **F-C2a, F-C2b, F-C2c pre-registered**, k = 1 each; all three REFUTED. The 72-cell grid IS priced as an order statistic (`ST.best_of_k_within`: observed -0.6466 against an across-target null -1.2522, **194% accounted**; `ST.split_half_transfer`: -0.1486, 23% of the ORACLE gain) and is never quoted without that price. The remaining 150 cells are descriptive map characterisation and spend no alpha. |
| 18 | 2026-09-21 00:37 | V | S31-V1 matched-K control for the set-matched readout ladder: 2 fold-clustered paired contrasts (exhaustive 8128-pair min vs argmin128; matched-K 128-random-pair min vs argmin128) + a 14-point x 2-family order-statistic growth curve (descriptive) | 2 | CA point cloud, 126 targets, ORACLE / NOT DEPLOYABLE throughout; NO built-chain claim | **V1 pre-registered in `s31/AUDIT_V.md` s0 (commit bac1fe1c) before the number existed**, k = 1 for its own purposes. The bar FIRED: at matched K the pair family is +0.0771 A WORSE (1.11x MDE, 5/5 folds). The growth curve is a DESCRIPTIVE decomposition -- nothing is accepted or rejected on it -- and the matched-K arm reports the MEAN over 8 draws, never the maximum (contract rule 10). |
| 20 | 2026-09-21 01:05 | V | S31-V3 adversarial audit of AVG_SEP at 126/126: 1 fold-clustered paired contrast (chain) + 1 cloud cross-check + 1 reproduction check + a geometry-instrumentation census (descriptive) | 1 | BUILT CHAIN primary, CA point cloud cross-check; both arms from the SAME row of the SAME job (lane D standing rule 2) | Applies lane F OWN registered falsifier (PREREG_S31_F.md s10.3) verbatim; k = 1, and the family was already registered by lane F, so this row is an independent RECOMPUTATION and spends no additional alpha. VERDICT REFUTED: +0.4609 A WORSE, 2.34x MDE, fold CI [+0.378,+0.565], 5/5 folds, 41W/85L. |

> **LANE C APPENDED ROWS 15-17 AT 2026-09-21 00:25 and did NOT recompute the running totals below**,
> because the total line is a read-modify-write that has already collided once tonight. Lane C's
> contribution is **k = 13 + 76 + 77 = 166**, of which the alpha-spending families are the 72-cell
> index grid (already priced as an order statistic) and nothing else; the 8 pre-registered
> falsifiers F-C1a/b/c/d, F-C2a/b/c and F-C3a/b are k = 1 each for their own purposes. Whoever
> next edits the totals should add 166 to the registered total and 72 to the alpha-spending total.

**k = 332 comparisons registered** as of 2026-09-21 00:10
(126 + 126 + 2 + 28 + 28 + 18 + 2 + 2 = 332, recomputed — the first version of this line said
187, which was wrong by 145 and is corrected in place rather than silently; a multiplicity
register whose own total does not add up is worse than none).

**α-spending families only: k = 60 + lane B's 146 = 206** (rows 4, 5, 7, 8, 10, 11, 12; lane B
counts rows 10-12 as α-spending in full even though row 10 is a decomposition of someone else's
published statistic and row 11's three primaries were pre-registered -- the conservative choice is
the honest one here). *Original line, left standing:* **α-spending families only: k = 60** (rows 4, 5, 7, 8). Rows 1–3 are instrument
characterisation and row 6 is a descriptive decomposition: nothing was accepted or rejected on
them, so they spend no α, and that judgement is recorded here in the open rather than left to
whoever would benefit from it.

At **k = 60** the adjusted multiplier is **4.183**, i.e. **1.49× the single-comparison MDE**, for
any cell of an exploratory family. Pre-registered primaries keep 2.8016.

*Lanes: append your row when you emit the family. A lower bound in the hundreds is what S30
ended with and it is not a number anyone can correct with.*
| 7 | 2026-09-21 01:13 | E | E5 magnitude control family: 4-point shrink curve c*y through the endpoint | 8 | both bases, 126 targets | ORACLE magnitude reference curve; no effect claimed from any single point |
| 8 | 2026-09-21 01:13 | E | E5 direction controls: energy-matched split and isotropic split of y, x2 bases | 8 | both bases, 126 targets | the matched null for "is mu special"; **built before the endpoint arms ran, launched after seeing them** -- stated in S31-L22 rather than left to be inferred |
| 9 | 2026-09-21 01:14 | E | E7: the 15 paired contrasts S31-L22 quotes, persisted to an artefact (rule 17), x2 bases | 30 | both bases, 126 targets | no NEW claim -- these are the contrasts already reported, given an artefact path; 5 are pre-registered, 10 are their matched controls |
