# PREREG — S32 LANE Q (QUANTUM / CVaR-VQE REFORMULATION)

Registered **before the first number**. Branch `s26`. Endpoint = **mean built-chain Cα RMSD,
n = 126, `tuning126`; production 3.2105 Å**. Cloud (3.0483) and set mean (3.5507) are different
objects and are never differenced against it.

Statistical rule for every arm below: `s24/stats_lib.compare(v, base, folds=..., names=...)`,
LOWER IS BETTER. **< 0.7× MDE → NOT A RESULT. 0.7–1.0× → NOT MEASURED. ≥ 1.0× with fold CI
excluding zero and ≥ 4/5 folds → RESULT.** SE printed beside every mean.

---

## Q0 — FALSIFY S31's TARGET-INVARIANCE CLAIM

S31 §5.1 claims the deployed CVaR-VQE Hamiltonian `E = _zrank(pool["sc"][o])`
(`core/pipeline.py:863`) is the same vector on every target to within **0.0407 max-norm**, and that
the target therefore enters the answer only through the readout's `P` and `W`.

### Q0-H1 (structure) — the residual is EXACTLY the tie pattern and nothing else

**MECHANISM.** `order = argsort(sc)` (`:757`), `top = order[:want]` (`:759`), `o = top[:128]`
(`:863`) ⇒ `sc[o]` is non-decreasing ⇒ `rankdata` of a sorted vector is the ramp `1..128` **except**
where values tie, where the whole block takes its average rank.

**PREDICTION.** For every target, `E` equals `zrank(midrank(sc[o]))` **bit-for-bit**, and
`max|E − ramp| > 0` **iff** the target has ≥ 1 tie in its top-128.

**FALSIFIER.** Any target where `sc[o]` is not non-decreasing; or where `E ≠ ramp` with zero ties;
or where `E = ramp` with ≥ 1 tie; or `max|E − ramp| > 0.0407 + 1e-9`.

**SELF-TEST THAT CAN FAIL (contract rule 5).** S31's own check was run on random floats, which never
tie. Mine is run on the **real 126 pools**, and additionally on a synthetic vector **constructed to
tie**, and on a **deliberately unsorted** vector, which must produce a non-ramp `E`. If the
deliberately-tied and deliberately-unsorted cases do not break the ramp, the test is decoration and
is discarded.

### Q0-H2 (paths) — no other target-dependent quantity enters the OPTIMISATION `min_p F(p)`

**PREDICTION.** The complete argument list of the optimisation is `(E, α, T, n, layers, iters,
seed)`. Of these, `α, T = VQE_LFO[fold % 5]` is **fold**-dependent (two distinct cells), and
`n, layers, iters, seed` are constants in `Config`. `dim = min(128, len(top))` is 128 on 126/126 or
the stage raises. Therefore `p` is a function of `(tie pattern, fold)` only.

**FALSIFIER.** A path — pool size < 128, `cfg.m`, duplicate handling, `tie_break`, quantisation —
by which any other target-specific quantity reaches `F`. Enumerated by reading `:846-877` and
checked numerically for `len(top)` and `dim` on all 126.

### Q0-H3 (pricing) — the 0.0407 residual buys less than 0.7× MDE

**0.0407 is not zero and must be priced, not dismissed.** At the deployed `α = 1, T = 0.3`,
`p* ∝ exp(−E/0.3)`, so a 0.0407 shift in `E` is a **14 % weight ratio change** on a tied block —
per-element non-negligible. The question is what it is worth at the endpoint.

**ARMS** (all at the fold's deployed `(α,T)`; `p*` = the exact hinged-Gibbs minimiser, since S31
proved the circuit is strictly worse 126/126 and the circuit adds only optimisation noise):

| arm | `E` used | readout |
|---|---|---|
| `Q0_TIE` | the real `E` (carries the tie pattern) | convex `average_weighted` |
| `Q0_RAMP` | the exact ramp `zrank(1..128)` — the target-independent constant | convex `average_weighted` |
| `Q0_TIE_SEL` / `Q0_RAMP_SEL` | the same two | selection (`consensus_medoid`) |

**PREDICTION.** `Q0_TIE − Q0_RAMP` is **< 0.7× MDE** on the CA cloud, and the selected medoid is
identical on ≥ 120 / 126 targets.

**FALSIFIER — and this is the one that would make S31 wrong.** If `Q0_TIE − Q0_RAMP` reaches
**≥ 1.0× MDE** on the CA cloud, the tie pattern is a real target-dependent channel inside the
Hamiltonian, S31 §5.1's conclusion is overstated, and the arm is **carried to the built chain**
(contract rule 0: a cloud result is a diagnostic until projected) and reported as the sprint
headline.

**CONTROL.** `Q0_PERM`: the same tie-block structure applied at **randomly relocated** positions of
the ramp — matched in number of tied elements and block sizes, unmatched in *where* the ties are.
This is the control matched to the operator's own space (contract rule 6): if `Q0_TIE − Q0_RAMP` and
`Q0_PERM − Q0_RAMP` are the same size, the channel is the *existence* of ties, not *which*
candidates tie, and it carries no target information. **≥ 8 draws, reporting the draw mean and the
draw-to-draw sd (contract rule 10), never the best draw.**

**DEPLOYMENT CONDITION.** None. Q0 is an audit of a prior claim; no arm here is a deployment
candidate.

---

## Q1/Q2 — A GENUINELY TARGET-DEPENDENT FORMULATION

Charter §16 requires all fifteen items **before compute**. Item 13, *classical equivalent*, is the
one that killed S31 and is decided first for every objective; an objective a classical algorithm
solves exactly and cheaply is reported as closed and not measured.

### Q1-T1 (THEOREM, registered as a derivation with a numerical falsifier)

**CLAIM.** For any weights with `Σ_x w_x = 1`, the readout objective `‖Σ_x w_x W_x − t‖²` depends on
the native `t` **only through the linear functional `g_x = ⟨W_x, t⟩`**, and therefore only through
`P_aff{W} t`, the projection of the native onto the affine hull of the candidate set. Consequently
the per-candidate quality vector `a_x = ‖W_x − t‖²` and the pool common mode `μ = X̄ − t` are the
**same object up to a known, native-free affine bijection**, and the sufficient statistic for the
entire readout problem has dimension `rank(aff{W_x}) ≈ 33`, not 128.

**WHY IT MATTERS.** S31 §20.3 states as *"the sharpest question this sprint produces, and it is
stated as open"* whether `â` (which candidate is better) and `μ` (how they are all wrong together)
are one requirement or two. If the claim holds they are one.

**NUMERICAL FALSIFIER (must be able to fail).**
1. Reconstruct `a` from `t` by the affine map and compare to the directly computed `a`: **must agree
   to < 1e-9 relative**, on all 126 targets.
2. Recover `P_aff t` from `a` alone by least squares and re-form `a`: **must agree to < 1e-9**. If
   the recovery is rank-deficient in a way that loses readout-relevant information, the claim fails.
3. Solve the simplex QP `min_w ⟨w,a⟩ − ½w'Bw` with the **true** `a` and with the `a` reconstructed
   from `P_aff t`: the two minimisers must give **identical** objective values to < 1e-9.
4. The identity `‖Σ w W − t‖² = ⟨w,a⟩ − ½w'Bw` is re-verified independently (S31: 1.66e-11).

**FALSIFIED IF** any of 1–4 fails, or if `rank(aff{W_x})` is not the dimension that controls the
recovery.

### Q1-A (readout posterior / CVaR over quality uncertainty) — classical-equivalent test first

**PREDICTION.** Under a Gaussian posterior `a ~ N(â, Σ)`, `CVaR_α` of the readout error over that
posterior is `⟨w,â⟩ − ½w'Bw + φ(α)·√(w'Σw)`, which is **convex on the simplex** (the first two terms
are the convex quadratic `‖Σw W − t̂‖²` by Q1-T1; the third is a norm). Therefore it is a **second-
order cone program with a unique global optimum, solvable classically in milliseconds**, and by
charter §16.13 it is **closed, not measured** — *unless* a cardinality constraint is added.

**FALSIFIER.** A demonstration that the objective is non-convex on the simplex for a legitimate
posterior, or that the posterior cannot be taken Gaussian without losing the mechanism.

### Q1-C (sparse convex combinations, `s`-of-`K`) — the only lead predicted to survive

**PREDICTION.** Choosing `s` of `K` plus weights is **cardinality-constrained** convex QP = best-
subset selection, which is NP-hard in general and is **not** a total ordering, so S31 §5.2's
obstructions are checked one at a time rather than assumed: (1) per-shot eigenvalue — a subset
bitstring `x` has the definite value `E(x) = min_{supp(w) ⊆ x} ‖Σ w W − t̂‖²`, so **diagonal in the
subset basis**; (2) not mean-field — `E(x)` is a genuine function of the bitstring; (3) dimension —
the register is over **subsets**, `2^K`, not over candidate index, `K`.

**PREDICTED KILL, registered in advance so it is not a post-hoc excuse.** `E(x)` needs `â`, and S31
§20.2 measured every native-free `â` at in-band `ρ` of **zero or the wrong sign**. By Q1-T1 an `â`
with real in-band skill **is** a structure estimate. So the expected verdict is *the combinatorial
problem is genuine and its objective is unknowable*, which is a different closure from S31's.

**FALSIFIER of that kill.** A native-free `â` with in-band `ρ > 0` at n = 126, or a formulation of
`E(x)` that does not require `â` at all.

### Q1-B (reconstruction branches) — owned by lane R; Q prices only the decision class

**PREDICTION.** If the branch decision is **per-target with a handful of branches**, it is
exhaustively enumerable and classically closed by charter §16.13. It escapes only if the branch is a
**per-residue** binary decision, giving `2^{n_res}` with a diagonal per-shot energy. Q reports which
of the two it is, from lane R's artefacts, and does not re-measure lane R's Ångströms.

### Q4 — is the five-bit structure an artefact, and can the information live in amplitudes?

**PREDICTION from Q1-T1.** The sufficient statistic is `r ≈ 33` **real numbers**, so a bit count is
the wrong currency; the five-bit result prices a *selection* alphabet, not the decision. Registered
measurement: the **ORACLE / NOT DEPLOYABLE** curve of emitted error against the number of retained
principal components `r` of `P_aff t` in the candidate basis, `r = 1 … 33`. This is a *pricing*
curve in real numbers, reported as ORACLE on every occurrence and never as a deployable arm.

**FALSIFIER.** If the curve is flat until `r ≈ 33`, the statistic is not compressible and the
"amplitudes can carry it" framing is dead. If it saturates at small `r`, the compressible direction
is named, with its ORACLE label attached.

---

## Multiplicity

Every comparison emitted by this lane is appended to `s32/MULTIPLICITY.md` **as it is emitted**,
tagged REGISTERED (the arms above) or EXPLORATORY (anything else). Registered arms: `Q0_TIE`,
`Q0_RAMP`, `Q0_TIE_SEL`, `Q0_RAMP_SEL`, `Q0_PERM`(≥8 draws), `Q1-T1`(1–4), `Q1-A`, `Q1-C`, `Q4`.

## What this lane will NOT do

Optimiser, depth, iteration or ansatz tuning (charter §36) — no objective in this lane has endpoint
evidence yet, so none of it is in scope.
