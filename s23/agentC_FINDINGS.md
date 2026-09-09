# WORKSTREAM C -- Sprint 23 findings (CVaR-VQE selector)

Pre-registered in `s23/PREREG_C.md` before either experiment read a single RMSD; the six-axis
fork list was sent to the coordinator before running, per Rule 0. Both experiments are complete.
No AMBER/OpenMM anywhere in this workstream, no benchmark inspection.

**One-sentence pillar statement, as BRIEF demands.** H_C1 is a genuine CVaR-VQE experiment: what
is quantum is the trained amplitude PATTERN that sets the intra-tail weights (an object no prior
closure covers); what is not quantum is the tail's MEMBERSHIP, which s22 D8 already proved is
always the classical top-m set. H_C2 does not run a VQE at all -- it conditions the classical
m-ladder directly, which is licensed only because s22 D8 proved that ladder IS what the CVaR
tail's membership always reduces to; the routing threshold itself is a plain classical scalar
over a native-free signal.

---

## 0. RESULT IN ONE LINE

**H_C2 is a clean, well-powered NULL, exactly as pre-registered pessimistically.** **H_C1 is NOT
a null -- it is a well-powered NEGATIVE result, and my own registered expectation ("trained ~=
untrained ~= classical bar, a clean null") is PARTLY WRONG.** Probability-weighting the CVaR tail
never beats the size-matched classical average at any tested temperature (0.0/0.1/0.2/0.5); at
the pre-registered PRIMARY (T=0.5) it ties the bar (null), but at both pre-registered SECONDARY
cells (T=0.1, T=0.2) it is SIGNIFICANTLY WORSE, and the harm is worst at intermediate
concentration -- a coherent, mechanistically explicable shape, not noise. **This is the "clean
and publishable negative" the BRIEF invited, but sharper than the form anticipated**: not
"training does nothing" but "the one door s22's theorem left open is closed empirically, and the
closure has a shape." Neither experiment moves the incumbent; both extend prior closures (s22 D8
for H_C1; s22 L7/L10 for H_C2) to new, previously-untested objects.

---

## 1. H_C1 -- PROBABILITY-WEIGHTED CVaR-VQE READOUT

### 1.1 Setup, exactly as pre-registered

Candidate-identity register (`s22/qcand_lib.Encoding`, reused as-is), `H|i> = E_i|i>` diagonal on
the shipped distogram Bayes-risk score, K<=500 -> n_qubits=9 (dim=512), exact `StatevectorCircuit`
via `core.quantum.run_cvar_vqe` (layers=3, iters=80, lr=0.15, restarts=1). alpha=0.15. T
PRIMARY=0.5, secondary T in {0.1, 0.2}. 4 seeds/target (0-3), n=126 dev targets, point-cloud
basis throughout, all 5 pinned folds.

Readout: `X_out = sum_{i in tail} (mass_i / sum mass) * X_i`, mass from `cvar_from_probs`,
against (a) size-matched classical top-m average at the trained tail's realised m [THE BAR],
(b) unweighted average over the identical tail SET, (c) the untrained circuit (`iters=0`,
identical seed) weighted the same way, at its own realised tail.

### 1.2 GATE 1 (soundness), reported before any RMSD

**100.0% pass, both arms, every T, all 126 targets x 4 seeds (2016 cells).** The trained-tail
candidate SET is element-for-element identical to `argsort(scores)[:m_realised]` on every single
cell, for both the trained and the untrained circuit. s22 D8's theorem reproduces exactly, as
required before any RMSD was read. No padding state ever entered a tail at alpha=0.15 (`pad_trained`
= `pad_untrained` = 0 everywhere) -- the sentinel penalty energy works as designed.

### 1.3 PRE-REGISTERED numbers (T=0.5 PRIMARY, T=0.1/0.2 secondary)

n=126, 4 seeds/target, point-cloud RMSD. `W` = probability-weighted average, `C` = size-matched
classical (uniform) average over the identical SET -- by Gate 1 this also equals arm (b), the
unweighted average over the trained tail, so (a) and (b) are reported as one column throughout.

| T | mean m (tr / un) | entropy bits (tr / un, max 9) | W_trained | C_trained (bar) | **PRIMARY diff** | verdict |
|---|---|---|---|---|---|---|
| **0.5 (PRIMARY)** | 66.0 / 78.1 | 8.85 / 5.81 | 3.0621 | 3.0408 | **+0.0213** SE 0.0117 MDE 0.0327, `\|eff\|/MDE`=0.65, CI_iid [+0.002,+0.048] excl.0, CI_fold [-0.001,+0.055] incl.0 | **NULL** (fold CI touches zero) |
| 0.2 (secondary) | 39.7 / 78.1 | 8.38 / 5.81 | 3.2285 | 3.1514 | **+0.0772** SE 0.0200 MDE 0.0560, `\|eff\|/MDE`=1.38, CI_iid [+0.039,+0.118], CI_fold [+0.038,+0.121], both excl.0 | **MEASURED -- SIGNIFICANT HARM** |
| 0.1 (secondary) | 29.4 / 78.1 | 8.16 / 5.81 | 3.3142 | 3.2068 | **+0.1074** SE 0.0232 MDE 0.0650, `\|eff\|/MDE`=1.65, CI_iid [+0.064,+0.155], CI_fold [+0.080,+0.142], both excl.0 | **MEASURED -- SIGNIFICANT HARM** |

**At every T, the probability-weighted average is at best tied with (T=0.5) and at worst
significantly WORSE than (T=0.1, T=0.2) the size-matched classical bar. It is NEVER significantly
BETTER, anywhere.** W/L at T=0.1 is 44W/80L; at T=0.2, 49W/77L; at T=0.5, 54W/72L -- weighting
loses on a majority of targets at every T.

**Control (c), the untrained circuit's own weighting at its own tail -- identical across T since
`iters=0` never sees T** (mean m=78.1, entropy 5.81 bits): `W_untrained`=3.2081 vs its own bar
`C_untrained`=3.0565, **diff +0.1517** [CI_fold +0.109,+0.211], MDE 0.1133, `|eff|/MDE`=1.34 --
**SIGNIFICANT HARM, larger than any trained cell.** So (c) does NOT tie (a): the untrained
control is significantly worse than its bar, more so than every trained condition. Per the
pre-registration's own stated reading rule, this is the OPPOSITE of "weighting is doing the work
and training is not" -- training measurably reduces the harm relative to an untrained circuit's
raw amplitudes at every T (0.152 -> 0.107 -> 0.077 -> 0.021 as T rises from untrained through
0.1, 0.2, 0.5), but never converts it into a net benefit.

**MECHANISM decomposition (training's effect on the WEIGHTED readout vs the UNIFORM readout at
the realised set, same T):**

| T | `W_tr - W_un` (training's effect on the weighted avg) | `C_tr - C_un` (m-realisation alone, uniform weights) |
|---|---|---|
| 0.5 | **-0.1460** [-0.232,-0.092] MEASURED, training helps | -0.0157 [-0.073,+0.015] NULL |
| 0.2 | +0.0204 [-0.070,+0.129] NULL | +0.0949 [+0.050,+0.149] MEASURED, harm |
| 0.1 | +0.1061, CI_fold excl.0 but CI_iid touches zero -> NULL | +0.1503 [+0.089,+0.219] MEASURED, harm |

At T=0.5 (the least-collapsed, most-uniform pre-registered cell), training SIGNIFICANTLY rescues
the weighted readout relative to the untrained circuit (-0.146), and the pure size-realisation
channel (uniform weights, trained-vs-untrained m) is null -- so the rescue is genuinely about the
WEIGHTS, not the set size. At T=0.1/0.2, that rescue is not established (both NULL or borderline),
and a meaningful part of `C_tr - C_un`'s harm is just the trained tail being SMALLER (m=29-40 vs
78) landing on a slightly worse part of the still-fairly-flat m-ladder.

**The effect is concentrated, not spread evenly (median vs mean gap, per project practice).**

| T | mean diff | median diff | top-10-targets' share of the total | worst target |
|---|---|---|---|---|
| 0.5 | +0.0213 | +0.0012 | **107%** (rest of the panel roughly cancels) | 2NDN +1.24 |
| 0.2 | +0.0772 | +0.0037 | 63% | 7QZV +0.93 |
| 0.1 | +0.1074 | +0.0064 | 55% | 9BAF +1.38 |

The mean is set by a small number of targets where weighting fails badly; the median effect is an
order of magnitude smaller than the mean at every T. This is the project's recurring
median-vs-mean signature and is reported rather than smoothed over -- the harm is real on
average but it is not "every target gets slightly worse," it is "most targets are close to
neutral and a minority fail hard."

### 1.4 EXPLORATORY ADDENDUM, dated 2026-09-08, NOT pre-registered

`s23/c1_explore_t0.py`, T=0.0 (pure CVaR, no entropy term), same alpha/seeds/architecture,
appended AFTER the registered primary/secondary numbers above were already on disk. Motivation
(disclosed to the coordinator before running): entropy stayed pinned near the 9-qubit maximum
(8.16-8.85 of 9 bits) across all three registered T's, so none of them is a strong test of
weighting in a genuinely non-uniform regime; T=0 forces real collapse (s22 A2's own finding).

    mean m: trained 1.7 (near-total collapse)  entropy trained 4.69 bits
    W_trained=3.4298   C_trained(bar)=3.3913
    PRIMARY diff = +0.0386  SE 0.0107  MDE 0.0299  |eff|/MDE=1.29
    CI_iid [+0.018,+0.060]  CI_fold [+0.018,+0.054]  both excl. 0   VERDICT: MEASURED -- HARM
    W_tr - W_un = +0.2217 [+0.181,+0.271] MEASURED -- training makes it WORSE here (sign flips
                                                        from T=0.5's rescue)
    C_tr - C_un = +0.3348 [+0.291,+0.386] MEASURED -- reproduces the already-closed s22 A2/D5
                                                        fact that an unregularised CVaR collapse
                                                        loses to a wide classical average, ~0.33 A,
                                                        close to (a touch above) the established
                                                        averaging-operator price (-0.28 to -0.31 A
                                                        independent of the energy, s21)

**This completes a coherent, hump-shaped picture across all four T's tested (0.0 exploratory,
then 0.1, 0.2, 0.5 pre-registered):**

    T           0.0     0.1     0.2     0.5
    mean m      1.7     29.4    39.7    66.0
    harm (Å)   +0.039  +0.107  +0.077  +0.021   (all >= 0; only 0.5 is a null, none is a win)

**The harm is smallest at both extremes and largest at intermediate concentration.** At T=0 there
is almost nothing to weight (mean m=1.7, often literally m=1, where weighted and classical are
IDENTICAL by construction) so weighting has little room to hurt. At T=0.5 the trained distribution
is itself nearly uniform (entropy 8.85/9 bits) so the weights barely differ from uniform. The
danger zone is the MIDDLE: a distribution concentrated enough to be genuinely non-uniform but not
concentrated enough to collapse to a near-point estimate is exactly the regime where deviating
from a uniform mean adds variance without adding any RMSD-relevant direction -- because (Gate 1)
the trained amplitude pattern, however non-uniform, is a function of the SAME energies the
classical rank order already uses, and nothing in this project's closure history gives it a
channel to carry information the rank order does not already have.

### 1.5 Disposition

**The pre-registered falsifier's "clean negative" framing anticipated (c) tying (a) -- i.e.,
training doing nothing. That is not what happened.** Training does something measurable and
mostly protective (it reduces the harm relative to an untrained circuit at every T, and at T=0.5
specifically rescues the weighted readout back to parity with classical) -- but it never converts
weighting into a net win, and at two of three pre-registered cells plus the exploratory T=0, the
weighted readout is SIGNIFICANTLY WORSE than the classical bar. **The one door s22 D8's theorem
left open is now closed empirically, with a mechanism (variance without direction, worst at
intermediate concentration) rather than left as "untried."** This is reported as damaging to my
own registered expectation (which predicted a flat null) precisely because BRIEF asked for that
kind of honesty, and it is a stronger, more informative result than the clean null I predicted --
a real, well-powered, mechanistically coherent negative, not an absence of signal.

---

## 2. H_C2 -- rg_z-CONDITIONED AGGREGATION WIDTH

### 2.1 Setup, exactly as pre-registered

Classical top-m ladder from `s21/results/poolgap.json` (avg_500/150/75/20/5/1), bit-verified:
`avg_75` mean = 3.0483380938795324, exactly the incumbent. `rg_z` extended to the full 126-target
dev instrument (`s23/qc_lib.build_rg_table`, `s21/rgsign.py`'s exact recipe -- verified
bit-identical, max|diff|=0.0, on the 75-target overlap with `s21/results/rgsign.json`). Router:
ONE fitted threshold `tau` on rg_z, TWO pre-declared rungs never fit (m_lo=75 for rg_z<tau,
m_hi=150 for rg_z>=tau, direction fixed a priori), nested leave-one-fold-out CV over the 5 pinned
folds.

### 2.2 PRIMARY result

    fixed m=75 (incumbent)             3.0483
    routed(rg_z), nested 5-fold CV     3.0496   (36/126 routed to m=150)

    routed - fixed:  +0.0012  SE 0.0196  MDE 0.0550  |eff|/MDE 0.02
    CI iid  [-0.0373, +0.0401]     CI fold [-0.0202, +0.0334]     W/L 18/18 (90 exact ties)
    VERDICT: NULL / NOT MEASURED

**The falsifier does not fire.** The effect (+0.0012 A) sits at 2% of its own MDE -- not merely
"not significant" but essentially exactly zero, on the full 126-target instrument with no
subsampling. This is as clean a null as this project's statistics can produce: SE is tight
(0.0196), so the null is a genuine absence of signal, not an underpowered non-result.

### 2.3 Secondary, reported regardless of sign

    rg_gap in place of rg_z:      diff -0.0017   |eff|/MDE 0.03   CI_fold [-0.0216,+0.0295]   NULL
    REVERSED direction (adversarial control): diff +0.0000, exactly.

The reversed-direction control is itself informative: nested CV, faced with a direction the
mechanism does not support, self-disables (fits `tau` outside the data range on every fold, i.e.
"route nobody") rather than manufacturing a false signal. That is the correct, conservative
behaviour for this router construction and rules out one way this result could have been an
artefact of the CV procedure rather than the feature.

### 2.4 Disposition

**NOT DEMONSTRATED**, joining s22 L7's four prior router failures as a FIFTH independent
construction, and consistent with s22 L10's own bound (a single global threshold's own
generalisation gap, 0.39 A at n~100/fold, is comparable to the entire 0.482 A routing ceiling).
rg_z is mechanistically different from every feature L7 tried (a compactness-DISAGREEMENT signal,
not a functional of the objective's own score distribution) and had never been tried as an
aggregation-width gate specifically -- so this null is new information, not a repeat, but it does
not move the needle. **rg_z's demonstrated skill (L27/L28: partial rho 0.34-0.38 against the
in-band ORDERING outcome) does not transfer to the AGGREGATION-WIDTH lever at this sample size.**
That is consistent with, and extends, s22's own diagnosis: the signal is real, the router class
the sample size supports is too weak to spend it, regardless of which native-free feature is used.

---

## 3. Artefacts and completion

All under `s23/results/`, atomic writes (`.tmp` + `os.replace`), config-derived content
(`config` block records alpha/T/seeds/layers/iters/lr/n_qubits inline), full-key-set completion
flags (`complete` requires `len(rows) == 126`, not a row count alone):

| file | role | complete | n |
|---|---|---|---|
| `rg_table126.json` | rg_z family, full 126-target extension of `s21/rgsign.py` | True | 126 |
| `c2_rgcond.json` | H_C2 primary + secondary + context | True | 126 |
| `c1_probweight_raw.json` | H_C1 pre-registered raw (T=0.5/0.2/0.1 x 4 seeds x trained/untrained) | True | 126 |
| `c1_probweight_analysis.json` | H_C1 pre-registered aggregated stats | True | -- (derived) |
| `c1_explore_t0_raw.json` | EXPLORATORY T=0.0 raw, dated 2026-09-08 | True | 126 |
| `c1_explore_t0_analysis.json` | EXPLORATORY T=0.0 aggregated stats | True | -- (derived) |

Code: `s23/qc_lib.py` (shared stats + rg_z), `s23/c1_probweight.py` + `c1_analyze.py`
(pre-registered H_C1), `s23/c1_explore_t0.py` + `c1_explore_t0_analyze.py` (dated addendum),
`s23/c2_rgcond.py` (H_C2, self-contained). `s23/results/` also holds other lanes' concurrent
artefacts (`agentA*`, `agentB*`, `d_scale*`, `gscale.json`) -- untouched, listed here only to note
they are not this workstream's output.

## 4. Scope limits, stated plainly

- H_C1's register (n_qubits=9, K<=500) is Sprint 22 Workstream A's exact setup, chosen for direct
  comparability with its closed findings (D1-D8) -- not the shipped pipeline's own
  `quantum_stage` (n_qubits=7, dim=128) or its `average_weighted` operator (full-support
  weighting, a different, already-existing object).
- alpha=0.15 was chosen to size-match the incumbent's own m=75 rung; a different alpha would
  realise a different tail size and was not explored as a systematic axis (declared, not hidden).
- H_C2 conditions the classical ladder rather than retraining the VQE's own alpha per threshold;
  this is licensed by s22 D8's set-equality theorem (the tail membership at any alpha IS the
  classical top-m set) and was declared as the FUNCTIONAL fork before running, not discovered
  as a shortcut afterward.
- Neither experiment touches AMBER or Legacy; both are native-free at inference (rr/nat_ca used
  only for the evaluation label, never for training or routing).
- **Note for whoever reads `core/pipeline.py` next.** `quantum_stage`/`average_weighted`
  (pipeline.py:818-883) already implement a DIFFERENT, wider probability-weighted average in
  production: full-support p_theta weighting over the whole filtered set (not CVaR-tail
  restricted), evaluated in `bench_results/fourcomponent_tuning126_w8.json`'s
  `rmsd_q_avg`/`rmsd_q_synth` vs `rmsd_u_synth` (uniform-weight ablation) -- also a null there
  (`quantum_synthesis_vs_uniform`: -0.013 [-0.081,+0.054]). That object was NOT re-tested here;
  H_C1 is specifically the CVaR-tail-restricted weighting s22 D8b named as untried. The two
  results are consistent in direction (weighting does not beat uniform) but are not the same
  experiment and should not be cited interchangeably.
