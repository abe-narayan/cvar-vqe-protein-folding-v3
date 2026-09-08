# Sprint 14 figures

Each figure answers one scientific question. Regenerate them all with:

    python -m s14.figures

Every figure is built from the sprint's own result JSONs in `s14/results/` — nothing is
hand-entered — and a figure whose input is missing is skipped rather than faked. Arms that
read a native quantity are labelled ORACLE in the figure itself, because the most dangerous
failure mode in this repository is an oracle diagnostic quietly becoming a headline.

| figure | the question it answers |
|---|---|
| `fig01_native_free_ladder.png` | **What does every native-free method actually emit?** Nothing beats the retrieval pipeline that consumes the very same windows — all six torsion arms lose, and every confidence interval excludes zero. Every arm passes the leakage guard, which replaces the native trace with noise and demands bit-identical output. |
| `fig02_error_coherence.png` | **At matched angular accuracy, does the SHAPE of the error matter?** Yes, by up to 0.93 Å at σ = 12°, and the σ needed to reach 2.0 Å ranges from 10.6° to 20.2° on coherence alone. But all ten real emitters sit in the near-i.i.d. band, so the Sprint 12 restraint surface stands and is mildly conservative — and the coordinator's hypothesis that real predictors have coherent errors is refuted in the figure's own caption. |
| `fig03_positional_cost.png` | **Where along the chain does a torsion error cost anything?** A symmetric hump peaked at mid-chain: the middle 40% carries 62.8% of the cost, the outer 40% carries 15.7%, and the terminal tenth is ~14× cheaper. Geometric, not empirical — a torsion at position p hinges two segments and the cost follows the lever-arm product p(n−p). This explains the Sprint 13 terminal-dropout correction mechanically and means coverage must be position-weighted. |
| **`fig04_MONEY_budget_saturation.png`** | **Does searching a GOOD objective harder produce better structures?** No. A 2,000-fold increase in search buys 0.17 Å, all of it by evaluation 300, while the objective improves monotonically and the best available structure improves from 2.907 to 1.532 Å. The selection gap grows to 2.04 Å. Read the three columns against each other: separately, each one misleads. |
| `fig05_objective_quality.png` | **Do any of these objectives rank structures where a search actually lives?** The structural objective reaches in-decile +0.161 under a uniform proposal against Legacy's +0.043 and raw AMBER's −0.088. The caption carries the retraction: a first reading of +0.370 was withdrawn because the proposal was drawn partly from the prior it scores, and a sequence-blind twin control is reported beside it. |
| `fig06_averaging_space.png` | **Does the averaging SPACE matter, and what does physical validity cost?** Coordinate averaging beats torsion averaging of identical windows by +1.024 Å, and projection onto ideal geometry **costs** 0.157 Å (18W/108L). The caveat is in the caption: the raw coordinate average is not a valid backbone, so that number is a price tag on physical validity, not a free win — and every torsion method pays it by construction. |
| `fig07_set_consensus.png` | **Should a well-ordered objective pick one structure or select a set to average?** A set. Consensus improves from 3.608 to 3.314 Å while the set members never improve — all of it error cancellation in coordinate space. It still loses to the incumbent by +0.110. This is also the classical control any quantum claim must beat, since a CVaR tail measured B times and averaged is exactly this operator. |
| **`fig08_MONEY_causality_decomposition.png`** | **What does each stage contribute, and how much room is left for VQE?** Aggregation is worth 3.4× more than the objective (−0.587 Å against −0.171 Å). A VQE's only possible contribution is to prepare a better distribution than the prior — the job the objective does — and that job is worth 0.171 Å, already banked classically. Two independent measurements say no optimiser can add to it. |
| `fig09_certified_optimum.png` | **Is the objective's optimum somewhere you would want to go, and which term puts it there?** Yes, finally: the structural objective's certified global optimum is 0.885 Å *better* than a random draw where the coarse-grained energy's is 0.139 Å *worse*. Sprint 13's "the optimum is in the wrong place" is a property of the energies, not the problem. The right panel is the inversion: adding distogram weight improves global ordering monotonically while making the certified optimum monotonically worse. |
| **`fig10_MONEY_transfer_collapse.png`** | **Is in-band discrimination unlearnable, or merely untransferable?** Untransferable, decisively. The overfitting gap is −0.0005 [−0.0011, +0.0000] — indistinguishable from zero — while the transfer gap is −0.3859 [−0.4366, −0.3320], 770 times larger. A *linear* pair potential already saturates the within-target problem at 0.986 against the 0.638 that 2.0 Å requires, so no amount of capacity helps. Pairs are matched on separation, so the null is 0.505, not the 0.525 that applies to unbinned tail statistics. |

## Design notes

Palette is the three-slot categorical set validated in Sprint 13 against this light surface
(blue `#2a78d6`, orange `#eb6834`, aqua `#1baf7a`): all-pairs CVD separation and
normal-vision floors pass. Aqua sits below 3:1 contrast against the surface, so every aqua
mark carries a visible direct label — the relief rule. Series identity is never colour-alone:
each figure carries a legend, direct labels, or both. Grids are recessive and drawn below the
marks. ORACLE arms are drawn in neutral grey so they can never be mistaken for a method.
