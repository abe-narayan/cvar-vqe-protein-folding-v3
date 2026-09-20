# s29_P_FINDINGS (Sprint 29, lane P: the projection price)

Standing format (S12 to S28): DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN,
every number with its artefact path, then "what damaged my own expectations" and "what I did not
do and why". Pre-registration `s29/PREREG_S29_P.md` (written 2026-09-20 00:10, before any
number). Code `s29/s29_P_scale.py`, tests `tests/test_s29_P.py`, results `s29/results/s29_P_*`.
Ledger entries `## S29-L<n> -- ... (date time, P)` in `s29/LEDGER.md`.

The lane's object: production's point cloud is 3.0483 A and its built chain is 3.2126 A, so the
last stage of the pipeline COSTS +0.1643 A (measured here on all 126, not quoted from the
record). That stage is the only one no S28 lane touched. The cloud is geometrically inconsistent
with the space it is projected into. The lane asks whether removing that inconsistency before
the projection changes the built chain, and what the ceiling of the whole rescale family is.

---

## DEMONSTRATED

P1. **The production point cloud reproduces bit-exactly in this lane's process on 126/126.**
`s29/results/s29_P_factors.json` (`n_cloud_mismatch` 0): every target's `rmsd_cloud` recomputed
through `s29_P_scale.py :: production_cloud` (the `s27/s28_B_prodcheck.py :: project_production`
recipe) differs from `s27/results/chain_rows.jsonl :: DIS` by **exactly 0.0**. Job
`s26/jobs_done/s29P_factors.json` exit 0, 20.0 s, peak RSS 0.319 GB. Consequence: every contrast
in this lane has both sides on one code path and one input, so the S28-L18 / S28-L27b
cross-path floor does not enter; the *input-perturbation* floor is measured separately (arm
FLOOR).

P2. **S23 L1's geometry reproduces on the current production cloud, to four decimals.**
Mean adjacent CA-CA distance **2.9614 A** (S23 L1: 2.961), radius of gyration **6.2061 A**
(S23 L1: 6.206), against the ideal covalent 3.80 and a native Rg near 6.60. The cloud is 22%
short in the bond and 6% short in the envelope: **the contraction is not a uniform shrink**, and
that is a fact about the current pipeline, not only about the S23 one.

P3. **The bond-consistency factor g is far from a global constant.**
g = 3.80 / (mean adjacent CA-CA distance), native-free: **mean 1.352, sd 0.318, range 0.998 to
2.141**, deciles 1.002 / 1.026 / 1.342 / 1.611 / 1.758 (`s29_P_factors.json :: per_target`).
A third of the panel wants essentially no rescale and a tenth wants more than 1.75x. This is
why the matched-random control (a derangement of the realised g across targets) is a real
control here and not a near-copy of BOND.

P4. **The posterior's implied envelope is larger than the cloud's, but only by 10%.**
s_SPAN (the factor that makes the cloud's Rg equal the Rg implied by the distogram posterior's
L1-Bayes median map) has **mean 1.102, sd 0.084**; s_ISO (the least-squares slope of the cloud's
map onto the posterior median map) **1.077, sd 0.073**; the restricted variant s_SPANr 1.099.
**s_SPAN and s_ISO agree at Spearman +0.989** but correlate with g only at +0.63 and +0.56:
the two native-free rescale signals are NOT the same quantity, and they differ by a factor of
about 3.5 in magnitude (1.10 vs 1.35).

*(sections below completed when the 126-target run lands)*

## ORACLE DIAGNOSTIC

## HYPOTHESIS

## REFUTED

## OPEN

## WHAT DAMAGED MY OWN EXPECTATIONS

## WHAT I DID NOT DO AND WHY
