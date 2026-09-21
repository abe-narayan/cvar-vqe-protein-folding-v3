# PREREG — S31 lane A (quantum / CVaR theory)

Written and committed **before any number in it exists**. Every clause below is prospective
unless it is explicitly marked `VERIFICATION` (a code-correctness check on a theorem, which
carries no evidential weight about the science) or `ORACLE` (diagnostic only, NOT DEPLOYABLE).

Author: lane A. Date of writing: 2026-09-20. Branch `s26`.

Endpoint basis, stated once and re-stated at every number:

* **CA point cloud** — production canonical **3.0483 Å**. All of lane A's own measurements are
  on this basis, because the derivation in §2 is an exact identity about a coordinate average
  and says nothing about the torsion projection.
* **built chain** — production canonical **3.2105 Å**. Lane A does not measure on this basis.
  A cloud-basis correction is translated with the project's measured **0.92** cloud→chain
  transfer coefficient, and never with the 1.16 of `operator-consumes-set-mean`, which maps
  set-mean to output and is a different map.

MDE = 2.8016 × SE, **per comparison**, with fold-clustered SE on the pinned 5 folds.
Below 0.7× MDE = **not a result**. 0.7–1.0× = **NOT MEASURED**.

---

## 0. The objects

`D = 2^7 = 128` the deployed quantum hypothesis set (`core/pipeline.py:837`, `cfg.vqe_qubits = 7`).
`W_x ∈ R^(n×3)` candidate CA coordinates, Kabsch-superposed onto the **pool medoid of the 128**,
a frame fixed once and exogenously. `t` the native in the same frame (ORACLE).
`e_x = W_x − t`, `a_x = ‖e_x‖²_F` (ORACLE). `B_xy = ‖W_x − W_y‖²_F` (**native-free**).
`E = _zrank(DIS)` the deployed diagonal energies. `w ∈ Δ^(D−1)` a readout weighting.
`(α, T)` the deployed leave-fold-out cell, `core/pipeline.py:113`
`VQE_LFO = {0:(1.0,0.3), 1:(0.25,0.3), 2:(0.25,0.3), 3:(1.0,0.3), 4:(1.0,0.3)}`.

---

## 1. A1 — the monotone-reweighting theorem (what the deployed state can say)

**Claim.** For `F(p) = CVaR_α(E;p) − T·H(p)` on the simplex with `T > 0`, the minimiser is unique
and is `p*_x ∝ exp( −(min(E_x, q) − q) / (αT) )` with `q` the self-consistent α-quantile. Hence
`p*` is a **non-increasing function of `E` alone**; every level set of `p*` is a prefix of the
`E`-order, and the state's entire structural content is the pair `(α, T)` applied to one sort.
At `α = 1` (folds 0, 3, 4) this is exactly `p*_x ∝ exp(−E_x/T)`.

* **A1-v (VERIFICATION).** The closed form agrees with a numerical simplex minimiser of `F` to
  `< 1e-8` in max-abs probability on 200 random `E` at each deployed `(α,T)`.
* **A1-e (REGISTERED, empirical).** On all 126 targets, with the deployed per-fold `(α,T)`,
  `n=7, layers=3, iters=50, seed=0`:
  * **predicted:** mean `KL(p_θ ‖ p*) < 0.10` bits, **and** the cloud RMSD of the `p_θ`-weighted
    average and the `p*`-weighted average differ by `< 0.02 Å` in mean.
  * **falsifier:** mean `KL > 0.50` bits **or** mean |ΔRMSD| `> 0.10 Å`. If falsified, the
    deployed stage's content is an **optimisation-path artefact**, not the objective, and that
    is the finding to report.
* **A1-d (REGISTERED, degeneracy).** At `T = 0` the argmin set of `CVaR_α(E;·)` is
  `{p : p_{x₀} ≥ α}` for a unique-argmin `E` — dimension `D−1`. Predicted consequence: the
  `T = 0` readout is seed-dependent with per-target cloud RMSD sd `> 0.15 Å` over 8 seeds,
  while the deployed `T = 0.3` readout has sd `< 0.05 Å`. Falsifier: `T = 0` sd `< 0.05 Å`.

## 2. A3 — the exact readout identity, and the Hamiltonian it forces

**Claim (identity).** For any `w ∈ Δ^(D−1)` and any fixed frame,

        ‖ Σ_x w_x W_x − t ‖²  =  ⟨w, a⟩  −  ½ · w' B w        (exact, no approximation)

with `a` ORACLE and `B` **native-free**. Equivalently `½ w'Bw = tr Σ_w`, the weighted dispersion.

Consequences registered as claims:

* the **native-free part of the exact readout objective is pairwise and enters with a MINUS
  sign** — at fixed candidate quality the readout should **maximise** weighted mutual spread;
* the corresponding mean-field Hamiltonian is `H[w] = diag(â) − B`, whose VMC local energy
  `E_loc(x) = â_x − (Bw)_x` is the **exact gradient** of the readout error. This is the
  physically meaningful non-diagonal `H` lane A proposes, and it is forced, not chosen;
* therefore `H = diag(zrank) − λ·W(similarity)` with `λ > 0` (attraction / consensus) carries
  the **wrong sign** relative to the derivation.

* **A3-i (VERIFICATION).** Max relative deviation of the identity `< 1e-10` over 126 targets ×
  200 random `w` each.
* **A3-frame (REGISTERED).** The deployed `block` (`Pt`) is **pairwise RMSD with per-pair
  optimal superposition**, which is *not* the fixed-frame `B` the identity requires. Predicted:
  substituting `block²·n` for `B` breaks the identity with median relative error `> 1%`.
  Falsifier: `< 0.1%` (in which case the distinction is immaterial and should be dropped).

## 3. The arms, and the one primary comparison

All arms emit a CA cloud over the **same** deployed top-128 set, scored by Kabsch CA-RMSD to
the native. Baselines:

| tag | arm |
|---|---|
| `PROD75` | uniform coordinate mean of the tie-safe DIS top-75 — **production**, cloud 3.0483 Å |
| `TOP128U` | uniform coordinate mean of the top-128 — matched-set control |
| `QVQE` | the deployed `p_θ`-weighted mean over the top-128 |
| `BOLTZ` | the `p*` closed form of A1 — the classical reduction of `QVQE` |

New arms:

| tag | arm | status |
|---|---|---|
| `QPORACLE` | `argmin_{w∈Δ} ⟨w,a⟩ − ½w'Bw`, convex QP | **ORACLE — NOT DEPLOYABLE** |
| `MEB` | `argmax_{w∈Δ} ½w'Bw` — the minimum-enclosing-ball centre; zero parameters, no score | native-free |
| `CAL` | `argmin_{w∈Δ} ⟨w,â⟩ − ½w'Bw` with `â` calibrated **leave-fold-out** by regressing `a` on `z(DIS)`; **γ = 1 exactly, no free parameter** | native-free at inference, LFO-trained |
| `GAM` | `argmin_{w∈Δ} ⟨w,z(DIS)⟩ − γ·½w'Bw`, `γ` on a fixed grid, chosen **leave-fold-out** | native-free at inference, LFO-tuned |

**PRIMARY REGISTERED COMPARISON:** `CAL − PROD75`, mean paired difference over n = 126,
cloud basis, fold-clustered SE, MDE = 2.8016×SE.
**Predicted:** between **−0.15 and +0.10 Å**; I do not predict the direction.
**Falsifier for lane A's Q2 direction:** `|Δ| < 0.7 × MDE` ⟹ **NOT A RESULT**, and the derived
non-diagonal Hamiltonian buys nothing at the readout on the cloud basis. I will report that
plainly rather than searching for a sub-arm that crosses.

Secondary, each with its own MDE:

* **A3-sign (REGISTERED, decisive for the brief's candidate).** The `GAM` grid is symmetric and
  contains negative γ (`γ < 0` = the brief's attractive/consensus `−λW`, `γ > 0` = the
  derivation's repulsive `−B`). **The derivation predicts the LFO-selected γ is > 0 on ≥ 4 of
  the 5 folds.** Falsifier: LFO-selected γ < 0 on ≥ 3 folds ⟹ the sign claim is refuted at the
  level of the deployed readout and I will say so.
* **A3-meb.** `MEB − PROD75`. **Predicted +0.3 to +1.5 Å (WORSE)**, because `a` is far from
  constant and the MEB support is the pool's extreme points
  (`consensus-is-outlier-avoidance`). If `MEB` instead **beats** `PROD75` by > 1.0×MDE the
  derivation delivers with no quality model at all, which would be the sprint's result.
* **A3-oracle (ORACLE — NOT DEPLOYABLE).** `QPORACLE` mean cloud RMSD. **Predicted < 1.2 Å.**
  This is the ceiling of the entire weighted-average readout class and prices charter §7C.

## 4. Controls

* **shuffled-B.** Every new arm re-run with `B` replaced by `P B P'` for a random permutation
  `P` (spectrum preserved, candidate correspondence destroyed). Any effect that survives is not
  the geometry. 8 draws, the **distribution** reported, never its maximum
  (`draw-controls-need-their-own-distribution`).
* **shuffled-score.** `z(DIS)` permuted within target, for `CAL`/`GAM`.
* **matched-null for γ.** `γ` drawn uniformly from the same grid rather than LFO-selected;
  the LFO arm is compared to that distribution, not to the grid's best cell.
* **zero-information control.** `â` constant (which is exactly `MEB`), a *plausible* control
  rather than a uniform one (`zero-information-control-must-be-plausible`).
* **seed control.** 8 circuit seeds for `QVQE` in A1-d.

## 5. What lane A will NOT do

* not reopen benchmark60, not regenerate folds or clusters;
* not select `γ`, `â`'s calibration, `α`, `T`, the frame, or the arm set on the 126-target
  native RMSD — `γ` and `â` are leave-fold-out on the pinned folds and nothing else is fitted;
* not report any per-target maximum over variants as a gain (`grid-oracles-are-order-statistics`);
* not describe any construction here as a quantum mechanism without discharging §11 of the
  charter explicitly. Lane A's working expectation, stated in advance, is that **every arm in
  §3 is classically reducible in O(D²)–O(D³)**, and the report will say so.

## 6. Declared prior

I expect the **theorems** to be the deliverable and the **arms to fail**: the sprint record says
consensus is capped at the pool's mode, dispersion terms measured +0.2059 Å worse (S30-L12), and
the transferable halfspace rule landed +0.472 Å worse than production (S30 THEORY §4.3). I am
registering the arms anyway because the derivation fixes the coefficient and the sign that all of
those experiments left free, and that is a different experiment. If they fail I will report a
closure, which the coordinator's brief states is the better outcome.
