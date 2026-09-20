# S29 THEORY, IN THREE PAGES

Consolidation for the report, written for a reader who has not opened `s29/THEORY.md`. No new
derivations. Every claim carries its ledger number; the derivation behind it is in the section
named beside it. ORACLE quantities (anything computed from a native) are labelled where they
appear and none of them chooses a deployable parameter.

---

## 1. THE BOUND (S29-L23, section 8)

Write any native-free operator's output as `A = c + u`, `c` production, `u` a displacement computed
from the data alone. With `e = t - c` the ORACLE error and `rho = cos(u, e)`,

    |A - t|^2 = |e|^2 (1 - 2 rho s + s^2),   s = |u|/|e|,
    minimised at s = rho:   **RMSD_achievable = RMSD_prod * sqrt(1 - rho^2)**.

Selection, re-weighting, gradient steps, basin averages, signed readouts, projections, quantum
tails: all are displacements, so the whole operator space collapses onto one number. Inverting it
at production's 3.2126 A built chain gives the price list:

    target      3.00     2.50     2.31     1.71
    rho needed  0.358    0.628    0.695    0.847        (a RANDOM shape field gives 0.140, S28-L23b)

Measured `rho` for every field the project has built: shipped objective's descent **-0.034**
(S28-L23b), every other S27 channel +0.034 to -0.006, the typicality axis **-0.058** (S29-L20),
the pool's first shape mode with a coin-flip sign (S29-L21). A field's sign accuracy `q` enters as
`|rho|(2q-1)`, so PC1's per-target `|rho| = 0.37-0.39` at `q = 0.52` is worth nothing.

    every field ever built here (|rho| <= 0.04)                 >= 3.210 A   (0.003 A gain)
    a RANDOM-strength field with the sign right every time      >= 3.181 A   (0.032 A)
    PC1 with a PERFECT ORACLE per-target sign                   >= 2.98  A   (0.23 A)
    the charter's target                                         2.50  A   needs rho = 0.628

**Assumptions.** (B1) the operator is native-free. **(B2) `rho_max <= 0.14` for every field
constructible from the present information -- LOAD-BEARING, and the only one a new source can
break.** (B3) rigid-body components are projected out and the identity is read to first order in
`s` (under 2% at `s <= 0.14`). (B4) `rho_max` is an expectation over the 126, so a field that is
excellent on some targets and reversed on others enters at `|rho|(2q-1)`.

**Falsifier, one measurement:** a native-free field whose mean cosine with `t - c` over 126 targets
exceeds 0.140 with the fold CI excluding it. That is meter number 2 and it runs in minutes. Lane D
is running the field survey now; the deciding statistic is the **signed mean over 126**, not
per-target magnitudes (0.24 to 0.53 per target is exactly what `|rho| = 0.37` with a coin-flip sign
looks like).

**What the bound says about the ceiling.** It is not the pool (its hull holds a 1.12 A point and
the shipped top-75's hull a 2.00 A point, lane O's ladder), not the circuit (a 27-parameter family
holds a 0.25 A structure on every target, S28-L26b), not the optimiser (it reaches the objective's
optimum on 126/126), not the readout. **It is the absence of a native-free vector with a cosine
above 0.14 to the direction that matters.**

---

## 2. THEOREM 2, AND WHAT SURVIVED THE ATTACK (S29-L7, amended S29-L29; section 2)

**The identity.** For any objective that sees the structure through its distance map,
`grad f(c) = Jc^T g` with `g` the per-pair coefficients, and exactly

    cos(-grad f(c), t - c) = -<g, r> / (||Jc^T g|| ||t - c||),   r = Jc (t - c).

An objective is locally informative **iff its per-pair force coefficients anti-correlate with
production's own signed per-pair error against the native**. Nothing about functional form,
normalisation, temperature, readout or ansatz enters.

**The theorem.** Decompose every pair quantity about the typical map: `d(c) = tau + a`,
`m = tau + b`, `d(t) = tau + n`. Under (A1) the marginal class, (A2) a proper risk, (A3) coherence
of `a` and `b`, and (A4) no channel correlating with `n`,

    E[<g, r>] = - sum w kappa var(a) (1 - beta),   beta = cov(a,b)/var(a),

**which contains no term in `n`**: the expected cosine of any marginal objective is second order in
the small quantities and carries information about the pool, never about the native.

| result | status |
|---|---|
| the exact identity | **stands** (lane D reproduced both sides independently, S29-L26) |
| the theorem's central claim (no term in `n`, hence second order) | **stands**; the bound rests on it |
| **corollary 2a** -- no function of the marginals, separable or not, can have an expected cosine whose size is set by `n` | **stands** |
| **corollary 2c** -- non-separability buys nothing; the part of `g` in `ker(Jc^T)` (45% of pair space at N = 12) is annihilated exactly | **stands** |
| **corollary 2b** -- the SIGN law, `E[cos] < 0` iff `beta > 1` | **WITHDRAWN at my own registered bar, both limbs (S29-L26, accepted S29-L29):** beta median 0.57 to 0.76 under four definitions, above 1 on 10 to 13%, sign agreement 48 to 50% inside a coin-toss CI of [41, 59] |
| the shrink warning (contract rule 20) | **stands and is now mechanical**: a positive cosine is purchasable with zero information by shrinking the target map toward typicality, so meter number 2 always prints its shrink signature (S29-L10) |

**The post-mortem (section 2.3b, S29-L29).** `beta < 1` is the pool's idiosyncratic 32% appearing
in `a` and not in `b`; substituted into the formula it predicts a *positive* cosine while the
measurement is negative, so a dropped term dominates. Two candidates: **(i) the linearisation** --
`phi' = 2F - 1` saturates, and at production `|2F-1| >= 0.5` on a median 53% of pairs (6 targets; the 126-target run is in flight), where the
coefficient is `w sign(a-b)` and the expectation stops being a function of `beta` at all; **(ii)
(A4) is an idealisation** -- in full `E[(a-b)(a-n)] = var(a) - cov(a,b) - cov(a,n) + cov(b,n)`, and
the measured sign reads `cov(a,n) > cov(b,n)`: **the pool's deviation tracks the native better than
the posterior's median map does.** The discriminating test (the `beta` law on unsaturated pairs
only, with the saturated subset as the matched control) is `s29/s29_T_beta_unsat.py`, built on lane
D's own implementation so that the author of the failed claim does not re-derive a favourable
variant. Reading (ii), if it wins, is the one pointer this sprint produced that none of its
closures touch, and it argues for operators that read **the pool's own dispersion** rather than the
posterior's marginals.

---

## 3. THE STABLE-RANK LAW AND ITS DESIGN RULE (S29-L11, section 3)

Exactly, `dF/dtheta_k = <psi|[A, Gt_k]|psi>` with `Gt_k` real antisymmetric and `Gt_k^2 = -I/4`;
for a 2-design state and a unit-spectral-norm observable,

    **Var[dF/dtheta] ~= tr(A_0^2)/D^2 = r_stable(A)/D^2**,   r_stable = ||A||_F^2 / ||A||_2^2.

The variance depends on `A` only through its Frobenius norm, which is basis-independent -- so a
diagonal matrix with the same spectrum decays identically, which is S28-L11's observation derived.
Three checks with no free parameter: S28's Gaussian graph to **7%** at n = 4 and n = 9; S28-B2's
unexplained "30 to 60x" for the kNN graph is exactly `M/k = 46x`; lane D's `J* = 85.7` against the
formula's 88.

**The design rule.** An off-diagonal term is gradient-visible at the deployed width iff its stable
rank grows with the register (parity at `n = 9` needs `r_stable ~ 8000`). No dense kernel can,
**centered or not** -- measured on 12 real pools, centering raises `lambda_2/lambda_1` from 0.138
to 0.465 but `r_stable` only from 1.04 to 1.59, and the decay gets **worse** (-2.30 per qubit
against -1.83); the signed agreement matrix is the same story. Any Gram of structural deviations is
capped by its rank at `3 N_res - 6 <= 42`. **The one operator class whose stable rank grows is a
LOCAL MIXER** (`r_stable = D/n`, slope exactly -1, `J* = 11.8`), which is meaningful only in a
configuration-space encoding, not in a candidate-index one (S29-L15, Q2).

**And the ground state of a centered Hamiltonian is the pool's principal contrast**, so any
probability readout averages its two poles back to the pool mean (section 3d). Lane O measured the
consequence: the ORACLE best global step along PC1 is exactly 0.000 A and the leave-fold-out step
is +0.0071 A worse than production (S29-L21).

---

## 4. THE Q1 REDUCTION, AND M6's CORRECTED STATUS (S29-L15, amended S29-L29; section Q1)

CVaR's derivative is `(E - q)/alpha` on the strict tail and **exactly zero above the VaR**, so at
the realised `m = 74` of `D = 512` the objective is flat on **437 of 511 simplex directions**. What
the entropy term selects inside that flat set is exact: by Sion's minimax theorem plus the
Rockafellar-Uryasev form of the LOWER tail,

    F* = max_t [ t - T log sum_x exp((t - E_x)_+/(alpha T)) ],   p*(x) ~ exp((t* - E_x)_+/(alpha T)),

uniform above the VaR and exponentially enhanced below it (the code asserts the `alpha = 1` limit
returns the Gibbs free energy). At the deployed cell: uniform `F = -4.5394`, `m = 92`; the exact
optimum `F = -4.7237`, `m = 29`; the trained circuit `F = -4.5610`, `m = 74.1`.

`E` is the standardised rank ladder on every target (S25 L17), so `p*` depends on `(alpha, T)` and
nothing else. **The deployed quantum stage is a target-independent weight profile over RANKS, and
its only endpoint channel is where the prefix cuts -- one number, `m`** (measured sd 6.74,
correlation +0.026 with chain length). That is the full mechanism for "the optimiser reduces the
objective on 126/126 and the structure does not move".

**M6, corrected.** My registered clause -- the deployed arm equals the fixed profile to within the
built-chain floor on at least 120 of 126 -- is **FALSE at the structure level** (exact set equality
on 3 to 10 targets, within the floor on 43, within 0.02 A on 77, worst 0.58 A apart) and **TRUE at
the endpoint** (-0.0097 A, 0.43x MDE). Lane D reproduced the profile arithmetic independently
(`t*` -1.5329 vs -1.534, `F*` -4.7221 vs -4.7237, `m*` 30 vs 29, PR 340.4 vs 342.3), so the
derivation is confirmed and only my clause was wrong. **The report's phrasing: the deployed stage
is not bit-reproducible by a fixed profile but is statistically indistinguishable from one** --
readout slack for the third time (S25 L15). M6 stands as an **endpoint-equivalence control**, which
is stronger than "a classical equivalent" because it needs no circuit, no optimiser and no
per-target computation at all.

---

## 5. THE TAIL-THEN-AGGREGATE LIFT (S29-L17 section 4.5; measured by lane B, S29-L25)

The deployed CVaR acts on a scalar per basis state and its tail is provably a prefix of the energy
order -- which is why the selection carries no information the classical sort does not. Lift it so
that the observable is **the tail's own coordinate average**, `R_alpha(p) = sum_x lambda_x(p) W_x`,
with the objective `f(R_alpha(p))`. Then `dR/dp_y = (W_y - W_boundary)/alpha`, the parameter-shift
chain costs the same `2P` evaluations, and the estimator needs a quantile plus a **mean** structure
rather than the full `2^n`-outcome distribution -- **more** device-realisable than the deployed
objective, not less.

**The set-equality theorem fails for it**, with a three-state, two-dimensional witness: the optimal
2-subset is the two *worst* states by every per-state criterion, because their errors cancel;
`V(S) = f(mean_{i in S} W_i)` is neither additive nor monotone. **This is the first formulation in
the project's history whose classical counterpart genuinely goes away.**

**Lane B measured it on real pools (S29-L25), exhaustively over all 124,750 pairs of the deployed
500-pool, and confirmed both clauses I registered:** the `f`-optimal pair is not the energy-order
prefix on **11 of 12** targets and not the two best by the per-state criterion on 11 of 12 either;
the `f`-optimal `m = 5` subset is non-prefix on **12 of 12** with a mean objective gap of +0.152.
**And the escape buys nothing:** those subsets' averages are worse than the prefix on 10 of 12 and
worse than production on 8 of 12 (n = 12, direction only). Exactly as registered: escaping the
theorem is **necessary and not sufficient**, and the order of work is to find an `f` that clears
meter numbers 3 and 4 *first*.

---

## 6. WHAT WOULD HAVE TO BE TRUE FOR THE BOUND TO FALL

Only (B2) can break, and only a channel with `cov(., n) != 0` given the distogram can break it.
Ranked by what the record says each carries:

1. **A better distance prior.** It moves `e` itself rather than `u`: -2.15 A per unit of prior
   improvement, concave, the only steep lever ever measured (S24 `priorladder`). Blocked by
   hardware and leakage, not by theory.
2. **The native's per-separation profile -- 11 to 14 numbers.** An ORACLE scorer that knows only it
   selects a top-75 that is **0.75 A better than production** (2.299 vs 3.048, `s12/obj_FINDINGS.md`
   section 4), while the deployable version of the same scorer is 3.163, i.e. worse. The
   lowest-dimensional named instance of what would break (B2); its leading component is compactness,
   the in-band axis, where the oracle correlation is +0.909 and achievable proxies reach 0.24 to
   0.37 (S29-L19, S29-L24). Cheapest falsifier written for lane M or P.
3. **A learned residual** with errors decorrelated from the predictor's: first order by
   construction, blocked by error coherence (+0.31 A at 0.688 sign accuracy when the mistakes are
   coherent, S19).
4. **A physics term on the emitted structure, in its FREE-energy form.** Outside the marginal class,
   so theorem 2 does not bound it; its single-point form is measured worse than random (S25 L16);
   the free-energy form is the only class with peptide-length precedent (S29-L1) and S8's `F_qh`
   machinery exists and was never finished.
5. **The pool's own dispersion**, promoted by reading (ii) of the post-mortem: it is a second moment
   and can predict the error's magnitude but not its sign, and both binding meter numbers need the
   sign -- unless `cov(a,n) > cov(b,n)` survives its test, in which case the pool, not the
   posterior, is the better estimate of the native's deviation and that is a first-order statement.
6. **A second, differently biased pool** (breaks (A3) only: second order; 31% independent only when
   0.76 A worse, S24 L2/L3), **an ESM-attention map** (redundant with a distogram that already
   consumes ESM-2 650M PCA-32; S17), **a joint over the same marginals** (corollary 2c: annihilated
   at the gradient by construction). None of these changes the order of the answer.

**In one line.** The sprint's theory says the ceiling is an information statement with a number on
it -- a cosine of 0.14 where 0.628 is needed -- and the two places to attack it are the prior's own
accuracy and the pool's dispersion, not the objective's functional form, not the Hamiltonian's
off-diagonal structure, and not the ansatz.
