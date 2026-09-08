# SPRINT 15 FIGURES

Rendered by `python -m s15.figures`, which reads the result JSONs in `s15/results/` and regenerates
every panel. Each entry below states what the figure shows, which artefact it reads, and what it must
not be read as saying.

**An n-guard applies, with three declared exceptions.** `_load()` refuses any result with `n < 100`, printing
`REFUSING <file>: n = k < 100. Smoke read, not a figure.` This exists because three figures in this
sprint were briefly rendered from 6- and 8-target smoke runs and looked exactly like finished results.
The project's history makes that dangerous rather than untidy: an 8-target read once put an arm at
2.792 Å where the full instrument gave 3.644 Å and **reversed the sign of the comparison**.

**The three exceptions are declared in code and each prints its own n on the panel**: fig02 (40
targets by design), fig05 (126 targets, but sharing a file whose file-level count is 40), and fig15
(9 targets / 27 paired cells). An earlier version of this file called the guard "hard"; it was
bypassed three ways until an adversarial audit found them, and this is the corrected statement.

**Conventions.** ORACLE arms are drawn hatched and labelled ORACLE in the tick label — a ceiling is
never allowed to look like a result. The incumbent (3.204 Å) is a labelled reference rule on every
accuracy panel. Value labels sit inside bars. Axes start at zero except where a caption says
otherwise (fig16 only).

---

| figure | shows | reads | n |
|---|---|---|---|
| **fig01_cascade** | the generative cascade's four stages G/S/A/F against the incumbent, with paired intervals | `cascade_combined.json` | 126 |
| **fig02_phase_diagram** | **error shape beats error magnitude** — three error models on one comparable axis (effective RMS), with the real distogram marked | `distacc.json` | 40 *(declared exception; prints its n)* |
| **fig03_feasibility** | Family B's feasible set excludes the native at every usable ε | `feasible.json` | 126 |
| **fig04_native_percentile** | where each objective ranks the truth, and whether it orders the pool at all | `feasible.json` | 126 |
| **fig05_bias_profile** | the distogram's systematic over-prediction, growing with sequence separation | `distacc.json` | 126 *(declared exception: shares a file whose file-level count is 40)* |
| **fig06_ceiling** | perfect distance knowledge is worth 0.611 Å, not the ≈1.95 Å the project assumed | `distgeo.json` | 126 |
| **fig07_robust_losses** | *withheld* — the run was stopped as superseded; the hypothesis was refuted at n = 126 in `align_fit` | — | — |
| **fig08_augment** | mixing solved conformers into the retrieval pool — every weight is harmful, and cross-validation chooses w = 0 | `augment.json` | 126 |
| **fig09_expand** | undoing the averaging contraction before projection — every native-free scale reference fails | `expand.json` | 126 |
| **fig10_errstruct** | *withheld* — surrogate destruction; the run was deferred for capacity and the guard refuses its n = 10 partial | `errstruct.json` | 10 |
| **fig11_realizability** | **the errors describe a consistent wrong structure** — the real prediction is 2.4× closer to realizable than matched noise | `coherence.json` | 126 |
| **fig12_fusion_law** | ~~a parameter-free law predicting fusion gain from native-free structural disagreement~~ **SUPERSEDED 2026-09-06 (Sprint 16 RETRACT).** The two-member Krogh-Vedelsby (1995) ambiguity decomposition. The figure plots `law_prediction`, which uses the ARITHMETIC mean of `r` where the identity needs the QUADRATIC mean; the scatter's offset from the diagonal is 54% that error and 46% a frame convention, not a prediction error. **Do not publish as drawn.** See `s16/retract_FINDINGS.md` | `coherence.json` | 126 |
| **fig13_scale** | *withheld* — the run was stopped; its ORACLE ceiling (~0.1 Å) had already closed the direction | — | — |
| **fig14_channels** | no available channel reaches the accuracy 2.0 Å requires (0.566 against 0.638) | `relayed.json` | 126 |
| **fig15_cvar_trade** | concentration buys the set mean and pays the set best; α is a diversity dial | `relayed.json` | 9 targets / 27 cells *(declared exception; prints its n)* |
| **fig16_gaps** | where the ångströms go — selection and projection lose ground, aggregation recovers part of it | `cascade_combined.json` | 126 |

---

## Two figures need their provenance stated

**fig14 and fig15 read `relayed.json`**, which holds tables **relayed from agent workstream findings
files rather than recomputed by the coordinator**. The file carries a `_provenance` field and a
`_source` line per block naming the originating findings file. They are drawn with a footer saying so.
Every other figure reads a JSON written directly by the module that produced it.

## What these figures do not show

- **No figure shows a method that beats the baseline**, because none does. fig01 and fig16 show the
  central negative; fig09 shows a native-free intervention failing with its ORACLE ceiling beside it.
- **fig02's x-axis is the *effective* RMS of the injected error**, not each model's nominal σ. Plotting
  against σ would confound shape with magnitude and destroy the panel's only point. The `sep_scaled`
  model's moment factor is √E[sep²]/E[sep] = 1.1331 and the `outliers` model's is 1.375.
- **fig06's 0.611 Å is an ORACLE ceiling** and is drawn beside the predicted arm at 3.644 Å for
  exactly that reason. It also carries a start-draw uncertainty of sd 0.132 Å — quote it as ≈0.6 Å.
- **fig11 and fig12's non-ORACLE arms are native-free**; the null ladder in fig11 is ORACLE by
  construction, since destroying a property of the error requires knowing the error.
- **fig16 is the one panel with a non-zero-based y-axis**, stated in the panel itself, because it
  shows differences between stages that a zero-based axis would compress into invisibility.

## Regenerating

```
python -m s15.figures          # all panels; skips or refuses anything not ready
```

Individual panels can be rendered directly, e.g. `python -c "from s15 import figures as F; F.fig_cascade()"`.
