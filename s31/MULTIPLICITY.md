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
| 1 | 2026-09-21 00:0x | D | D-B projection determinism: per-target pass-1 vs pass-2 bit-identity | 126 | built chain, 126 targets | determinism assertion only; no effect claimed |
| 2 | 2026-09-21 00:0x | D | D-B 1e-14 input perturbation vs unperturbed reprojection | 126 | built chain, 126 targets | instrument property (noise floor); no effect claimed |
| 3 | 2026-09-21 00:0x | D | D-B reprojection vs s29 canonical rows, and vs s27 DIS rows | 2 | built chain, 126 targets | agreement check; no effect claimed |

*(Rows 1–3 are instrument characterisation, not hypothesis tests: they carry no α budget because
nothing was accepted or rejected on them. They are registered anyway, because the rule is that
comparisons are registered when emitted and the decision about whether they spend α is made
here, in the open, rather than by whoever would benefit from the answer.)*
| 4 | 2026-09-20 23:57 | E | E1 direction: cos/corr(mu_hat, mu) levels + 2 matched controls, 2 pair masks | 28 | within-target cosine, 126 (106 long-range) targets | **1 pre-registered primary** (uncentered cos, all pairs, vs the permutation control) named in `PREREG_S31_E.md` §2 before any number; the other 27 are its robustness surface and carry no separate claim |
| 5 | 2026-09-20 23:57 | E | E2/E3 applied correction arms: 14 arms x {built chain, CA cloud} vs PROD | 28 | both bases, 126 targets | **5 pre-registered** (E2/E3 native-free; E2o/E3o/class-ceiling under AMENDMENT 1) + 3 mandatory matched controls; the ORACLE decomposition arms are bounds, not hypothesis tests |
| 6 | 2026-09-20 23:57 | E | E3 orthogonal SSE decomposition of y against mu and mu_hat | 18 | ORACLE energy/SSE shares, 126 targets | descriptive mechanism decomposition; **no effect accepted or rejected on it** -- it is the explanation of an effect S30 already measured |
