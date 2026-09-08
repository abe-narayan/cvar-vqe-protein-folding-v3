# Sprint 13 figures

Seven figures, each answering one question. Regenerate them all with:

    python -m s13.figures

Every figure is built from the sprint's own result JSONs in `s13/results/` — nothing is
hand-entered — and a figure whose input is missing is skipped rather than faked. Arms that
read a native quantity are labelled ORACLE in the figure itself, because the most dangerous
failure mode in this repository is an oracle diagnostic quietly becoming a headline.

| figure | the question it answers |
|---|---|
| **`fig01_money_gradient_prediction.png`** | **Does a molecular energy's Pauli-weight spectrum predict its variational gradient variance?** Yes, with no free parameter. The left panel is predicted against measured over 104 cells per model; the right panel is the ratio on a linear scale, and it is the one to read — the left axis spans ~48 decades because raw AMBER's dynamic range does, and a log axis that wide flatters any fit. Median ratios: Legacy 0.997, conditioned AMBER 1.001, raw AMBER 0.880 and broad. |
| **`fig02_money_accuracy_ladder.png`** | **What does every native-free method actually emit, against what the space contains?** Nothing reaches the 1.594 Å the space holds; simulated annealing on Legacy is worse than random sampling. The right panel is the control that matters: the torsion prior's apparent advantage is a helix artefact, a zero-information constant α-helix reproduces it to within 0.10 Å, and on non-helical targets the prior's advantage vanishes. |
| `fig03_representation_ceiling.png` | **How much structure does the discrete torsion space hold per qubit spent?** At ~24 live qubits it contains a 1.594 Å answer, against a 3.213 Å retrieval pipeline. The marked X is the audit's correction: from a *random* descent start the identical search reaches only 1.982 Å, so 0.388 Å of the ceiling is the privileged oracle start. ORACLE diagnostic. |
| `fig04_objective_validity.png` | **Do the energies rank the native, and do they rank where a search actually lives?** No. Legacy's rank correlation inside the low-energy decile is +0.043. Raw AMBER spans 16 decades, and a *monotone* log compression leaves ρ unchanged while collapsing the range to 1.8 — the damage is in the ordering, not the scale. |
| `fig05_pauli_spectra.png` | **Is the all-atom force field a more global observable — and was the first answer real?** After conditioning, yes: AMBER peaks at Pauli weight 3 against Legacy's weight 1. The right panel is why the first answer was not: raw AMBER's top-10 configurations hold a median 99.6% of its variance, and a delta spike has Binomial(m, ½) weight for arithmetic reasons. |
| `fig06_budget_curve.png` | **Does optimising the objective harder produce better structures?** No. On nine fully enumerated targets Legacy improves to 300 evaluations then degrades, and its certified global optimum is worse than its own best intermediate. The ORACLE arm, searching on true RMSD, descends monotonically — that is what a good objective looks like. |
| `fig07_restraint_surface.png` | **How good would torsion information have to be, how good can we predict it, and where do the gaps fall?** The channel reaches 1.52 Å at σ = 12° with full coverage; the incumbent pipeline is equivalent to σ ≈ 29°, so a sequence-only builder must beat 29° merely to tie. The best leave-fold-out predictor achieves σ = 67.7° and 3.770 Å. The right panel carries a correction to the coordinator's own guidance: **terminal** dropout, which is where a chemical-shift predictor actually declines, is 0.40–0.50 Å *cheaper* than uniform dropout, so uniform understates the channel. |

## Design notes

Palette is the validated three-slot categorical set (blue `#2a78d6`, orange `#eb6834`, aqua
`#1baf7a`), which passes the all-pairs CVD and normal-vision separation floors on this light
surface. Aqua sits below 3:1 contrast against the surface, so every aqua mark carries a
visible direct label — the relief rule. Series identity is never colour-alone: each figure
carries a legend, direct labels, or both. Grids are recessive and drawn below the marks.
