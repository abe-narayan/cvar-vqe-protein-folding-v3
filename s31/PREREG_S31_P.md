# PREREG — S31 lane P: **substitute the closed-form `p*` for `run_cvar_vqe`'s `p` and score the endpoint**

**Written 2026-09-21, before the first number of this lane.** Lane P exists to run one
decisive experiment. Registered here: the arms, the primary comparisons, the falsifier in
Ångström on the built chain, the family size for `s31/MULTIPLICITY.md`, and my honest prior
odds — all before any measurement.

---

## 1. What is being tested, and why there is no third outcome

Lane L (`s31/LIT_L.md` §L1.1, ledger **S31-L4**) proved that the deployed CVaR free energy

```
F(p) = CVaR_alpha^low(E; p) - T * H(p)            core/quantum.py:993
```

is **convex in `p`** and has a **closed-form global minimiser** pinned by one scalar:

```
s*  = argmax_s { s - T * log sum_i exp( (s - E_i)_+ / (alpha*T) ) }     # 1-D, concave
p*_i  proportional to  exp( (s* - E_i)_+ / (alpha*T) )                  # hinged Gibbs
```

and that the shipped `run_cvar_vqe` is **strictly worse in 12/12 synthetic cells**, by an
**expressivity** gap (iterations buy nothing; depth and restarts do).

So the shipped CVaR-VQE is a lossy approximate solver for a convex program that has an
analytic solution. Substituting `p*` for the circuit's `p` therefore reads:

- **endpoint improves** → the quantum layer is a softmax plus a one-dimensional root-find,
  and the honest architecture is classical;
- **endpoint worsens** → the circuit's *inability* to optimise is the active ingredient,
  which is S20's law (`concentration-is-wrong-when-discrimination-binds`) arriving from the
  objective side, and is a real, publishable negative about the architecture.

Both are results. This prereg commits to reporting whichever fires.

## 2. `p*` is native-free — the trace, asserted before the run

`p*` is a function of `(E, alpha, T)` only.

| input | where it comes from | native? |
|---|---|---|
| `E` | `_zrank(pool["sc"][o])`, `core/pipeline.py:864` | **no** |
| `pool["sc"]` | `debias.score_risk(risk, GRID, D)`, `core/pipeline.py:731` | **no** |
| `risk` | `distogram_risk(seq, fold, ...)` — target sequence + the **held-out-fold** distogram model | **no** |
| `D` | pair distances of the retrieved pool windows | **no** |
| `o` | `order[:128]`, `order = argsort(sc)` | **no** |
| `alpha, T` | `VQE_LFO[fold]` = `{0:(1.0,0.3), 1:(0.25,0.3), 2:(0.25,0.3), 3:(1.0,0.3), 4:(1.0,0.3)}`, a **leave-fold-out** table fixed in S25 and not touched here | **no** |

The only place a native is read in the whole pipeline is `label()` (`core/pipeline.py:1140`),
which runs after every structure is final. `p*` is therefore **deployable**, exactly as the
VQE's `p` is. **No parameter of `p*` is tuned here on native RMSD, and `VQE_LFO` is used as
shipped.** This table is asserted in code (`s31/s31_P_substitute.py`, `--check-native-free`).

## 3. The arms — six, scored at the **built chain**

There are **three** readouts (`S31-L1` annotation), and `p` enters two of them.

| arm | stage | readout | operator |
|---|---|---|---|
| **A** | quantum OFF | — | `coordinate_average(W[top75])` → `project` — **the canonical 3.2105 Å** |
| **B** | ON, `p_vqe` | **selection** | `o[consensus_medoid(block, p_vqe)]` → `project` |
| **C** | ON, `p*` | **selection** | `o[consensus_medoid(block, p*)]` → `project` |
| **D** | ON, `p_vqe` | **convex** | `average_weighted(W[o], block, p_vqe)` → `project` |
| **E** | ON, `p*` | **convex** | `average_weighted(W[o], block, p*)` → `project` |
| **F** | ON, uniform | **convex** | `average_weighted(W[o], block, 1)` → `project` — the **shipped** zero-information ablation (`core/pipeline.py:1116`), used as built, not rebuilt |

`o = order[:128]`, `block = pairwise_rmsd(W[o])`, `alpha, T = VQE_LFO[fold]`,
`p_vqe = run_cvar_vqe(E, alpha, T, n=7, layers=3, iters=50, seed=0)` — the **deployed**
settings (`Config.vqe_layers=3, vqe_iters=50, vqe_seed=0, vqe_qubits=7`).

**The stage is OFF in production**, so a quantum arm must recover whatever the stage costs
merely to reach parity with the classical path it replaces. **See AMENDMENT 1 — the
+0.2260 Å figure I was briefed with is withdrawn and the size of that deficit is one of the
things this lane measures.**

---

## AMENDMENT 1 — the parity deficit is withdrawn, and two predictions are registered

**Filed 2026-09-21, still before the first number of this lane: no arm had been run when this
was written.** The coordinator traced the **+0.2260 Å** deficit I was briefed with (S31-L2,
"a quantum arm must recover 0.226 Å to reach parity") and it is **measured on the wrong
operator**: it is `circ_opt` = `circ_l1_i80` from `s27/s28_A_amp.py`, the **affine amplitude**
readout (signed weights, leaves the convex hull), not the shipped `average_weighted` convex
readout that arms D/E/F use. **The figure is withdrawn from this prereg.** This is the third
readout conflation of the sprint and the second a lane has caught, which is itself the reason
this lane scores all three readouts rather than one.

**What replaces it, and what is still missing.** Lane A measured the shipped operator this
hour: the quantum synthesis costs **+0.0178 Å on the CA point cloud** against production
(3.0661 vs 3.0483) and **+0.0125 Å** against the matched uniform-128 control (3.0536). **The
built-chain figure does not exist.** It is exactly `S3_D_minus_A` in §5 below.

> **Registered consequence: `S3_D_minus_A` (arm D − arm A) is promoted from a secondary
> context row to a REPORTED HEADLINE ROW in its own right.** It is the honest built-chain
> cost of the shipped quantum stage — the component the charter names as the spine — on the
> basis the endpoint is actually reported on, and this project does not currently have it. It
> is *not* the `p*` question and will not be presented as one. Family size is unchanged at
> k = 9; nothing is added, one row is relabelled.

**Lane A's prediction, registered here so that I am scored against it rather than able to
claim it afterwards.** Lane A measures `KL(p_theta || p*) = 0.930 bits, TV = 0.378`, yet its
two readouts differ by **−0.0044 Å at 0.088× MDE (NOT A RESULT)** on the cloud while differing
**per-target by 0.102 Å in absolute value**. *The optimisation gap is real and cancels in the
mean.* Lane A therefore expects my substitution to land **near zero on the cloud**, and holds
that the open question is whether it also lands near zero on the **built chain** — which,
given the projection's ~1e13 amplification at branch-degenerate points (S31-L4), is genuinely
a different question.

> **This sharpens P1/P2's reporting obligation, and the obligation is registered now.** A
> paired mean is not a sufficient report of this experiment. For both primaries I will emit
> the **per-target |Δ endpoint| distribution** — mean, median, p90, max, and the count above
> the 0.0107 Å chain floor — beside the paired mean. If the mean is null but the spread is
> ~0.1 Å, **the honest finding is "the substitution changes the answer on most targets and the
> changes cancel", which is a different sentence from "`p*` makes no difference", and I will
> write the first one.** A null mean with a large spread will not be reported as "no effect".

**My prior, updated by the corrected baseline and stated before the run.** Item 4 of §7
(~85% that every quantum arm stays worse than arm A) was conditioned on a 0.226 Å deficit
that does not exist. With the true cloud deficit at +0.0178 Å I revise it to **~60% that arm
D is worse than arm A on the built chain**, and I note the direction is now genuinely open:
the stage is nearly free, so a small `p*` effect is no longer automatically swamped. §7 items
1–3 are unchanged; **item 3 (~2:1 that `p*` does not improve on the VQE) I explicitly do not
revise**, because lane A's cloud measurement supports it.

## 4. Endpoint basis, and the instrument's own floor

- **Basis: BUILT CHAIN.** Canonical production = **3.2105 Å**
  (`s29/results/s29_O_chain_rows.jsonl :: item=prod`). The CA point cloud is 3.0483 Å and is
  **a different object**; every number in this lane states its basis in the same sentence.
- **The projection is deterministic but chaotic** (S31-L4, defect D-B): no RNG, but the λ=0
  multi-start argmin is decided at a ~1e-7 objective spread among branches ~1e-1 apart,
  amplification ~1e13. **Therefore all six arms are projected in the SAME process, from the
  SAME stored pool cache, with arm A re-projected from the stored canonical cloud
  (`s29/results/s29_O_structs/<pdb>.npz :: prod`) rather than quoted.** Registered before the
  run: if arm A's re-projected mean does not reproduce 3.2105 Å to ≤ 0.0107 Å, the lane
  reports the re-projected value as its own basis and says so.
- **No built-chain claim below 0.0107 Å.**

## 5. Primary comparisons, and the falsifier in Ångström

Paired at target level, n = 126, `s24.stats_lib.compare`, fold-clustered CIs on the **pinned**
folds. `d = a − b`, negative = `a` better. **MDE = 2.8016 × SE, per comparison.**

- **P1 (primary): E − D** — convex readout, `p*` vs VQE `p`, built chain.
- **P2 (primary): C − B** — selection readout, `p*` vs VQE `p`, built chain.
- **P3 (descriptive, not a test): the per-target disagreement count** — how often
  `consensus_medoid(block, p*) != consensus_medoid(block, p_vqe)`. Lane L's caveat is that the
  selection readout is `argmin(P·p)` and therefore **piecewise-constant**, so a large TV
  distance in `p` need not move the medoid at all. **If the disagreement count is low, P2 is
  uninformative and I will say so rather than report a null**; the convex arms then carry the
  experiment.

### The falsifier, stated in Ångström on the built chain

> **Lane L's "improves" branch is FALSIFIED if `E − D ≥ 0` at ≥ 1.0× MDE** (i.e. `p*` is
> measurably *worse* than the circuit's `p` through the shipped convex readout).
> **Lane L's "worsens" branch is FALSIFIED if `E − D ≤ 0` at ≥ 1.0× MDE.**
> **If `|E − D| < 0.7× MDE`, the substitution is a NULL**; if `0.7–1.0× MDE`, **NOT MEASURED**.
> **And if `|E − D| < 0.0107 Å` the comparison is inside the instrument's reprojection noise
> and is reported as such regardless of what the MDE says.**

The same three-way gate applies to P2.

The correctness of `p*` itself is falsified, and the lane aborts, if **any** of:
`|F(p*) − φ(s*)| > 1e-6` (strong duality), `F(p*) > F(p_vqe)` on any target (global optimality),
or an independent mirror descent on the full 128-simplex reaching below `F(p*)` by > 1e-6.

### Secondary comparisons (context, each k = 1, none primary)

`B − A`, `C − A`, `D − A`, `E − A`, `F − A`, `E − F`, `D − F`.

## 6. Family size for the multiplicity register

**k = 9** — 2 primary + 7 secondary. Appended to `s31/MULTIPLICITY.md` when the family is
emitted, not at the end of the sprint.

## 7. Honest prior odds, registered before any number

The coordinator registered **~3:1 that `p*` is NOT MEASURED against the VQE arm at the
endpoint**, because the selection readout is piecewise-constant. Mine, and it differs by
readout because the two readouts are not the same experiment:

1. **P2 (selection): ~3:1 that it is NULL or NOT MEASURED.** I agree with the coordinator here,
   and for his reason. I expect the medoid to disagree on a **minority** of targets — my point
   prediction is **20–45 of 126** — and a comparison that is identically zero on two thirds of
   the sample has to carry a large effect on the rest to clear its own MDE.
2. **P1 (convex): ~60/40 that it IS measurable**, because `p` enters continuously there and
   L1.1 measured `H(p*) > H(p_vqe)` in 12/12 cells — `p*` is the *more entropic* distribution,
   so arm E should sit **closer to the uniform arm F** than arm D does. That is a directional
   prediction and I register it: **I expect `E` to lie between `D` and `F`.**
3. **Direction: ~2:1 that `p*` does NOT improve on the VQE** (i.e. `E − D ≥ 0`). Reason, and it
   is not a hunch: S31-L2's annotation records that the convex optimum over the deployed set,
   converged, is **3.0522 Å (CA cloud) against production's 3.0483** — the best convex
   reweighting the objective can find **is** the uniform average it already emits. If the
   objective's own convex optimum is not better than uniform, then solving the objective
   *better* has no reason to help, which is S20 pointing at the objective.
4. **~85% that every quantum arm (B–F) remains worse than arm A**, because the stage is OFF in
   production for a measured reason and needs +0.2260 Å just to reach parity.

**If (3) fires, the honest headline is the negative — the circuit's failure to optimise is the
active ingredient — and I will lead with it.** I am registering (3) precisely so that I cannot
later present it as the outcome I was hoping for.

## 8. Standards this lane is held to

- `benchmark60` is **not** reopened. Nothing here touches it.
- Folds and clusters are **read** (`ST.pinned_folds`), never recomputed.
- No deployable parameter is tuned on native RMSD. `VQE_LFO`, `lam=0.3`, `m=75`,
  `vqe_qubits=7`, `layers=3`, `iters=50`, `seed=0` are all used **as shipped**.
- Every ORACLE quantity is labelled ORACLE in the same sentence as its number. The only
  ORACLE read in this lane is the native CA trace used to score RMSD, after every structure
  is final.
- `x.get(k) or 0` is not used anywhere; missing keys raise.
