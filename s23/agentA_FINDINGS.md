# WORKSTREAM A — ENSEMBLE AGGREGATION — FINDINGS (Sprint 23)

`s23/PREREG_A.md` filed before any experiment ran. `s23/agentA_run.py` computed all three
experiments in one pass over the 126 dev targets (`s23/results/agentA.json`, complete, self-check:
`avg_75` reproduces the pinned incumbent 3.0483 exactly, and `avg_500/150/75/20` reproduce
`s21/results/poolgap.json` bit-for-bit). `s23/agentA_report.py` produced the stats package
(`s23/results/agentA_report.json`, `s23/results/agentA_report.txt`). Completion flag
`s23/results/agentA_COMPLETE` requires the full key set, not a row count.

**HEADLINE: ALL THREE PRE-REGISTERED FALSIFIERS FIRE. No arm in any of the three experiments beats
the incumbent `avg_75` (3.048 Å) with a CI excluding zero on the achievable side. This is a clean
negative result across the workstream's whole priority list, and it extends the project's now
six-times-repeated pattern (routing, hedging, diversity selection, weighted CVaR tails, and now
cluster-restriction and weighted averaging) of real oracle headroom that no declared native-free
rule can reach.**

---

## EXPERIMENT 1 — CLUSTER-RESTRICTED AVERAGING: FALSIFIER FIRES, BUT THE MECHANISM IS REAL

At every (m, k) tested — m ∈ {75, 150, 250}, k ∈ {2,3,4,5} at m=75, k=3 at m=150/250 — **every
achievable choice rule (`size`, `score`, `combined`, `random`) is at or worse than `avg_75`**, most
with CIs excluding zero on the WORSE side:

| cell | rule | mean diff vs avg_75 | MDE multiple | 95% CI (iid) |
|---|---|---|---|---|
| m75_k3 | combined | **+0.142** | 1.00x | [+0.042, +0.244] |
| m75_k3 | random | +0.439 | 1.99x | [+0.288, +0.598] |
| m150_k3 | combined | +0.092 | 0.65x (ns) | [-0.005, +0.190] |
| m250_k3 | combined | +0.124 | 0.71x (ns) | [+0.001, +0.254] |
| m75_k3 | **ORACLE** | **-0.369** | **2.74x** | **[-0.467, -0.282]** |
| m150_k3 | **ORACLE** | **-0.430** | **2.65x** | **[-0.547, -0.321]** |
| m250_k3 | **ORACLE** | **-0.484** | **2.20x** | **[-0.646, -0.340]** |

**The oracle ceiling is real and grows with m** (−0.26 to −0.48 Å, 2.2–3.5× its own MDE, CI
excluding zero at every cell, 84–109 wins out of 126) — confirming the basin-splitting mechanism the
brief predicted genuinely exists on this pool. **No achievable choice rule reaches it.** `size` is
the least bad (0.5–0.9× its own MDE, mostly not significant) but never beats the incumbent; `score`
and `random` are actively harmful and clearly significant.

**A damaging and counter-intuitive result, worth flagging on its own.** At m75_k3 the real
`combined` rule (which typically excludes ~3–4 genuine geometric outliers from the top-75, keeping a
71-member coherent cluster) is **worse than the incumbent by +0.142 Å (CI excluding zero)**. The
matched control `random_same_size_partition` — which removes the SAME NUMBER of members but chosen
at random rather than by genuine structural dissimilarity — is statistically flat (+0.005, 0.12× its
own MDE, dead null). **Removing real geometric outliers from an already score-filtered top-75 hurts
more than removing an equal number of random members.** Read together with `error-coherence-decides-
correctors` (memory) and S21 L7's diversity-selection wreckage, the likely mechanism is that those
few outlier members carry decorrelated (i.i.d.-like) error that CANCELS in the average, and a
genuine-structure clustering step selectively strips exactly that cancelling error out, concentrating
whatever systematic bias remains in the bulk cluster. **Caveat on the control's own power**: because
the excluded group is only 3–4 of 75 members, `random_same_size_partition`'s "largest random group"
is nearly `avg_75` itself by construction, so its near-zero result is expected on construction
grounds and is weaker evidence than its flat CI alone suggests — it rules out "any 71-of-75 subset
average looks like this," not "clustering per se is inert."

**Disposition: REFUTED as an achievable lever. The mechanism (oracle headroom) is real and is the
FIFTH instance of this project's positive-control/negative-result pattern.** Not promoted; the named
k-medoids fork was not needed because the negative result is already clean at the declared
hierarchical-clustering construction, and the oracle/achievable gap is the finding, not a borderline
call that would motivate a second clustering algorithm.

---

## EXPERIMENT 2 — WIDTH SWEEP UNDER CLUSTERING: THE OPTIMUM MOVES, BUT NEVER TO SOMETHING BETTER

| m | avg_m mean | avg_m worst | avg_m fail% | clust_m mean | clust_m worst | clust_m fail% | clust−avg |
|---|---|---|---|---|---|---|---|
| 20 | 3.091 | 7.88 | 13.5% | 3.273 | 8.05 | 16.7% | +0.182 (1.37×) |
| 75 | **3.048** | 8.07 | 12.7% | 3.190 | 8.45 | 13.5% | +0.142 (1.00×) |
| **150** | 3.072 | 7.02 | 14.3% | **3.140** (clust_m optimum) | 8.00 | 12.7% | +0.069 (0.42×, ns) |
| 500 | 3.396 | 6.91 | 18.3% | 3.300 | 8.49 | 15.1% | −0.096 (0.40×, ns) |

**Plain `avg_m`'s optimum stays at m=75 (bit-exact reproduction of S22 L5).** **`clust_m`'s optimum
moves to m=150 — the hypothesis's directional prediction is correct** — but **`clust_m` at its own
best m (3.140) is still worse than plain `avg_75` (3.048)**, and `clust_m` is worse than or
statistically tied with `avg_m` at its OWN matched m everywhere except m=500, where the difference is
a non-significant tie (0.40× MDE). **Worst-target and failure-rate columns do not favour clustering
either**: `clust_m`'s worst-target RMSD is equal to or higher than `avg_m`'s at 6 of 8 widths (up to
+1.0 Å worse at m=100/250), and failure rate is a mixed bag with no consistent direction.

**Disposition: the falsifier's literal text ("optimal m stays at 75") does not fire, but the
substantive conclusion it was testing for does — clustering does not produce a competitive width
anywhere on the ladder.** The shift to m=150 is a second-order fact about where clustering does
LEAST damage, not a route to improvement.

---

## EXPERIMENT 3 — WEIGHTED AVERAGING: TOTAL CLOSURE, EVEN AT THE IN-SAMPLE ORACLE

| family | in-sample-oracle param | in-sample vs avg_75 | nested-CV achieved vs avg_75 |
|---|---|---|---|
| score-softmax | c=2.0 | −0.0038 (0.12× MDE, ns) | **+0.0018 (0.09× MDE, ns)** |
| rank-power | **p=0.0 (= uniform = incumbent)** | 0.0000 exactly | **+0.0052 (0.62× MDE, ns)** |
| distance-to-medoid | **c=1e6 (≈ uniform = incumbent)** | 0.0000 exactly | **+0.0018 (0.27× MDE, ns)** |

**All three pre-registered falsifiers fire.** This is a stronger closure than Experiments 1–2: for
`rank` and `dist`, **the grid search over all 126 targets with full native information (the
leakage-ceiling in-sample oracle) selects the UNIFORM-WEIGHT endpoint of its own grid as the best
point** — there is no oracle headroom to be unreachable, because within these three declared weight
families **the optimum genuinely is uniform averaging.** `score`-weighting shows a hint of oracle
headroom (−0.004, but at 0.12× its own MDE — not even close to significant) that evaporates entirely
under honest nested CV (+0.0018). **Weighted averaging, in the three families this experiment
declared, is not merely unreached — it has nothing to reach.**

**Disposition: REFUTED, cleanly, at both the achievable and in-sample-oracle level.** Not promoted.

---

## VIRTUAL BOND LENGTH — THE CONTRACTION TRADE, REPORTED PER BRIEF §CRITICAL CONTEXT

Physical target 3.805 Å; incumbent `avg_75` already at 2.961 Å.

- **Experiment 1** cluster arms sit at 3.03–3.36 Å (LESS contracted than the incumbent, because
  restricting to one cluster removes some of the averaging-induced smoothing) — but they are still
  worse on RMSD, so this is not a case of "winning by contracting further"; if anything the opposite
  concern applies (less contraction, still worse), which argues against contraction being the whole
  story for this instrument's optimum.
- **Experiment 2**: `clust_m`'s virtual bond tracks `avg_m`'s closely at each m (both drop toward
  ~2.5–2.7 Å at m=500 as the pool widens and averaging deepens); no arm buys RMSD by contracting
  further than its plain counterpart at matched m.
- **Experiment 3** weighted arms were not separately vbond-instrumented (all ties at/near the uniform
  point by construction), so there is nothing to report beyond the incumbent's own 2.961 Å.

**No arm in this workstream wins by trading geometry for Cα-RMSD — every arm that changes RMSD in
either direction does so without a compensating vbond story, because none of them beat the
incumbent in the first place.**

---

## WHAT THIS MEANS FOR THE SPRINT

All three of Workstream A's priority experiments are **closed negatively**, joining a now
well-established family of "oracle headroom exists, no native-free/fixed-rule construction reaches
it" results. The one piece of positive information is methodological rather than an arm to ship:
**Experiment 1's oracle ceiling (0.26–0.48 Å, growing with m) is the largest CLEANLY MEASURED
headroom this workstream found**, and it is a genuinely new object (never measured before this
sprint per the `s21`/`s22` search). If a future lane wants to chase cluster-choice specifically, the
k-medoids fork named in `PREREG_A.md` and a feature set built from `spread(cluster)` /
`n_clusters`-style geometry (in the spirit of S22 B4's candidate-set geometry features, which were
"the least harmful family tried" for m-routing) would be the next attempt — but S22's finite-sample
bound (L10/C1) argues the same n=126 ceiling applies here too, so this should be expected to fail for
the same structural reason rather than a fixable engineering gap.

**Nothing from this workstream is promoted.** `avg_75` (3.048 Å, point cloud) remains the best
point-cloud aggregation arm found across Sprints 21–23.

---

## ARTEFACTS

- `s23/PREREG_A.md` — pre-registration (dated 2026-09-08, unedited)
- `s23/agentA_run.py` — driver, one pass per target, all three experiments
- `s23/agentA_report.py` — stats package (iid + fold CI, MDE, W/L, worst-target, nested CV)
- `s23/results/agentA.json` — raw per-target rows, `complete: true`, n=126
- `s23/results/agentA_report.json` / `agentA_report.txt` — full numeric report
- `s23/results/agentA_COMPLETE` — completion flag, full required-key set, self-check against the
  pinned incumbent constant (3.048) and against `s21/results/poolgap.json` (bit-exact on avg_500/
  150/75/20)
