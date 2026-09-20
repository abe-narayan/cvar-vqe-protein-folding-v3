# PRE-REGISTRATION -- S30 lane T (theory)

Registered **before** any array was loaded. Every number below is a prediction with a bar; the
verdict rule is stated per clause. Written 2026-09-20 (see the ledger entry's timestamp for the
`date` capture in the same command as the append).

Instrument: `s8/generate_univ/<pdb>.npz`, the 126 pinned universes. Fields used: `rr` (ORACLE
CA-RMSD of every length-n window to the native), `order` (stable argsort of -`sim`, the deployed
BLOSUM key), `W` (window CA coordinates), `nat_ca`. Nothing here chooses a deployable parameter;
every ORACLE quantity is labelled.

---

## A. WHY THESE THREE AND NOT OTHERS

The charter's lead L8 asks where 5.56 of 7 bits went. The question presumes the register is a
**channel**. My derivation (s30/THEORY.md sections 1-2) says it is a **codebook index**, and that
the value of a bit is set by the codebook's order statistic, not by the bit count. That
reformulation is only worth anything if the order statistic is measured, so these three
measurements are the ones that decide whether my sections 1-2 are arithmetic or rhetoric.

M3 is additionally the rank-collapse pre-check the coordinator required before any compute is spent
on a lifted objective: S29's first non-diagonal Hamiltonian failed because its spectrum was a
near-rank-one projector, and my Theorem T1b says the lifted (tail-then-aggregate) objective's
reachable tails live in a halfspace class whose dimension is the **stable rank of the candidate
feature matrix**, not its nominal dimension. If that stable rank is ~1, the lift is a relabelled
one-dimensional sort and must not be built.

---

## M1 -- THE ORDER-STATISTIC LADDER AND THE VALUE OF A BIT

For each of the 126 targets, within the deployed BLOSUM top-500 pool, compute the ORACLE
best-of-first-`N` by `rr` for `N = 1, 2, 4, ..., 512` (`N` = 2^R, R = 0..8.97), using the deployed
pool ORDER (so R = 7 is exactly the production top-128 the readout sees). Mean over 126.

**P1a (functional form).** The ladder is fitted by `D(R) = a + c * 2^(-R/gamma)` with
**R^2 >= 0.99** and **gamma in [2, 6]**.
*Falsified* if R^2 < 0.97 or gamma outside [1.5, 8].

**P1b (the marginal value of the 7th bit).** `-dD/dR` at R = 7 is in **[0.05, 0.14] A/bit**.
Rationale: S29-L44 measured the 7-bit best1_top128 arm at -1.0670 A total, an *average* of
0.152 A/bit, and a concave ladder must have a marginal below its average.
*Falsified* if the measured marginal at R = 7 is outside [0.03, 0.20].

**P1c (the floor is not the binding term at R = 7).** The fitted floor `a` satisfies
**a <= 1.6 A** and `D(7) - a >= 0.4 A`, i.e. at the production register width the pool is still
far from its own floor and bits are still worth something -- the diminishing-returns regime has
NOT been entered.
*Falsified* if `D(7) - a < 0.2 A` (bits already exhausted) or `a > 2.0` (floor above the charter's
target, which would close the candidate-index allocation outright).

## M2 -- WHAT THE 3,252 BITS OF RETRIEVAL CHOICE ACTUALLY BOUGHT

Matched control: a uniformly random 500-subset of the SAME universe (20 draws per target, seeded
per target). Compare the ORACLE best-of-500.

**Definition, registered now so it cannot be chosen after the fact.** Retrieval's *realised*
information is the **search-equivalent bits**
`b_ret = log2(N* / 500)` where `N*` is the number of uniformly random draws from the universe whose
ORACLE best-of-N matches the BLOSUM-500's. This is the currency the readout spends, so it is the
only currency in which the 3,252 bits of choice can be priced. Where the random ladder saturates
before matching, `N*` is read off the fitted random ladder and the extrapolation is declared.

**P2a.** `b_ret <= 3.0 bits`. The key's Spearman with true RMSD over the universe is +0.066
(`s8/generate.py` RESULT 2); a rank correlation that small cannot move a minimum order statistic
far.
*Falsified* if `b_ret > 4.0 bits`.

**P2b (the size of the gain in Angstroms, as a check on my Gaussian-copula reasoning).** The
BLOSUM-500's ORACLE best beats a random 500's by **0.10 to 0.30 A** on the mean over 126.
Derivation registered: under a Gaussian copula at rho = 0.066, conditioning on the top 2.9% of the
key shifts the normal-scored `rr` by `rho * E[Z | top 2.9%] = 0.066 * 2.27 = 0.150 sd`, and the
minimum of 500 shifts by the same 0.150 sd because the conditional sd is sqrt(1-rho^2) = 0.998.
With sd(`rr`) of order 1 to 2 A that is 0.15 to 0.30 A.
*Falsified* if outside [0.0, 0.45] A, and **note the sign**: a NEGATIVE value (the key worse than
random on the minimum) is a separate and reportable outcome, not a mere failure.

## M3 -- RANK COLLAPSE: THE PRE-CHECK BEFORE ANY LIFTED OBJECTIVE IS BUILT

Per target, superpose all 500 pool windows onto the pool medoid (the deployed
`coordinate_average` frame), centre, and form
- `A_coord` (500, 3n): the superposed coordinates, column-centred;
- `A_dist` (500, npairs): the pair-distance map with `min_sep = 2`, column-centred.

Report `r_stable = ||A||_F^2 / ||A||_2^2` (= the participation ratio of the covariance spectrum)
and `lambda_1 / sum(lambda)` for both, mean over 126.

**P3a.** `r_stable(A_coord)` is in **[1.5, 6]** and `r_stable(A_dist)` in **[1.2, 4]**. Rationale:
S29 section 3c measured Gram-of-deviations stable ranks of 1.04 (uncentred) to 1.59 (centred) on 12
pools, and `pool-error-is-68-percent-common-mode` says 68% of the squared error is a single shared
direction.
*Falsified*, and this is the clause that matters, if `r_stable(A_dist) > 10` -- in which case the
lifted objective has genuine multidimensional freedom and Theorem T1b's cap is loose.

**P3b (the decision rule, registered).** If `r_stable(A_dist) < 2.0`, then by T1b the
tail-then-aggregate lift's reachable-tail class is, to within one effective dimension, a **single
one-dimensional sort** and I will report the lift as CLOSED at the encoding level -- the same
verdict S29's near-rank-one Hamiltonian earned, reached before compute is spent rather than after.
If `2.0 <= r_stable(A_dist) < 10`, the lift is a `ceil(r_stable)+1`-parameter family and I will
report the exact reachable-set cap in bits against the 300.6 bits of set choice.

---

## B. WHAT I AM NOT CLAIMING

- Nothing here is an Angstrom claim about a deployable operator. M1's ladder is ORACLE throughout,
  by construction: it is the *codebook's* geometry, which is what the design question needs, and it
  is not reachable native-free (S29's bound, assumption B2).
- `b_ret` is a search-equivalence, not a mutual information. The two differ and I say which is which
  in THEORY.md section 1.3.
- M3 measures the pool's feature geometry, not any Hamiltonian's spectrum. It bounds the lift's
  reachable class; it does not price the lift in Angstroms.

---

# ADDENDUM 1 -- M4, THE QUADRIC CEILING (registered 2026-09-20, before the run)

Asked for by the coordinator after M3 closed the halfspace lift: *what is the ceiling of the
second-moment (quadric) class, before anyone builds it?* Same ceiling-before-compute discipline
that M3 just paid for.

**Construction.** Per target: superpose all 500 pool windows onto the pool medoid ONCE (an
EXOGENOUS frame -- declared, because a tail-dependent frame is an unstated operator). Feature map
`Psi` = the centred pair-distance map; project onto its top `k = 6` principal directions (`k90`
measured at 5.61 in M3b) and standardise -> `Z (500, 6)`. Fix `m = 75`. The emitted structure is
the plain coordinate mean of the selected 75 in the fixed frame; the score is ORACLE CA-RMSD to
`nat_ca`.

Four classes, and note their sizes at fixed `m`:

| class | order | free parameters | class size at fixed m |
|---|---|---|---|
| PREFIX | the deployed `DIS` energy | 0 | **exactly 1 set** |
| HALFSPACE | `<g, z_x>` | 6 | O(D^6) |
| QUADRIC | `z_x' G z_x + <g, z_x>` | 27 | O(D^27) |
| FREE | any 75-subset | -- | C(500,75) = 300.6 bits |

`HALFSPACE` is the `G = 0` slice of `QUADRIC`, so a deficit for QUADRIC at matched budget is a
statement about SEARCHABILITY, not about containment, and will be reported as such.

**Budget matching, because this is exactly the `grid-oracles-are-order-statistics` trap.** Each of
HALFSPACE, QUADRIC and FREE gets the SAME `K = 5,000` ORACLE-scored draws. The matched null is the
FREE arm: best-of-5,000 uniformly random 75-subsets. Any class that does not beat FREE at the same
K has bought nothing but an order statistic.

**P4a.** QUADRIC's ORACLE best is within **0.10 A** of HALFSPACE's at matched K.
*Falsified* if QUADRIC beats HALFSPACE by more than 0.20 A -- in which case the second moment is a
real and findable enrichment and lane Q should build it.

**P4b.** FREE (best-of-5,000 random 75-subsets) is within **0.15 A** of PREFIX. Rationale:
`operator-consumes-set-mean` (`d_out = 1.16*d_set_mean + 0.04*d_set_best`) plus the common-mode
identity (S23-L9: 68% of the pool's squared error is a shift every member shares, and a set mean
cannot touch it), so set means concentrate and 5,000 draws cannot separate them.
*Falsified* if the gap exceeds 0.30 A.

**P4c.** All three searched classes beat PREFIX by **less than 0.6 A**, i.e. far less than the
1.31 A the same pool yields to a single-member ORACLE argmin (M1's fitted floor). Rationale: the
mean is the wrong operator to spend selection on, which is the law's own content.
*Falsified* if any class beats PREFIX by more than 1.0 A.

**What M4 cannot say.** Every arm is ORACLE. M4 measures the CLASS's reachable ceiling, not any
deployable operator, and no native-free rule is proposed here that selects within any class.

---

# ADDENDUM 2 -- M5 AND P5, REGISTERED BEFORE LANE D REPORTS ITS COMBINATION

## P5 -- LANE D's FIELD-COMBINATION CEILING, PREDICTED FROM A DIFFERENT INSTRUMENT

For `k` fields with signed cosine `rho_0`, exchangeable pairwise correlation `c`, and the
`rho`-vector aligned with the Gram's leading eigenvector,
`rho_comb^2 = k rho_0^2 / (1 + (k-1)c)`, and `(1 + (k-1)c)/k = s_1`, the PC1 variance share of the
field Gram (the eigenvalues sum to `k`). Hence

    rho_comb = rho_0 / sqrt(s_1)          -- the COUNT cancels; only the Gram's rank enters.

At `rho_0 = 0.1128` (lane D's DISTPOT, S30-L5):

**P5a** the 21x21 field Gram has stable rank in **[1.3, 3.5]** and `s_1` in **[0.45, 0.80]**.
**P5b** the **leave-fold-out** ORACLE-optimal combination reaches `rho_comb` in **[0.11, 0.22]**,
point estimate **0.15**; it does NOT reach 0.30 and does NOT clear the 0.358 that 3.00 A needs.
*Falsified* if the LFO combination exceeds 0.25.
**P5c** the IN-SAMPLE all-21 figure will be inflated to **0.25-0.45** and must not be read as a
ceiling (21 coefficients on 126 targets; `C^{-1}` amplifies exactly the directions the fields do
not share).

The orthogonal estimate 0.33-0.36 corresponds to `s_1 = 1/11 = 0.09`, i.e. a Gram of full rank 11.

## M5 -- THE COMBINATION CEILING FROM THE POOL's GEOMETRY

> **DECLARED DEFECT, AND IT IS MINE. P5d BELOW WAS WRITTEN AFTER M5 HAD ALREADY RUN.** The heading
> on this block originally read "registered before the run" and that was false: I wrote
> `s30_T_combo.py`, ran it on all 126, saw the answer, and only then wrote this addendum. P5d
> therefore records what I *expected* and is **NOT a pre-registration**; it carries no evidential
> weight as a passed or failed prediction, and the fact that it comes out FALSIFIED must not be
> read as a registered falsification. P5a/P5b/P5c above ARE prospective -- they concern lane D's
> combination, which has not reported -- and only those three are to be scored.
>
> I am leaving the clause in rather than deleting it, because a deleted wrong prior is
> unauditable and because the direction of my error is itself information: I expected the span to
> bind and it does not.

Any native-free operator that reweights or selects pool members emits
`u = sum_x w_x (W_x - c)` with `sum w_x = 1`, so it lies in the span of the pool's deviations about
the production average `c`. Measure `||Pi_k e||^2 / ||e||^2` for the oracle error `e = t - c`, with
`Pi_k` the projector on the top-`k` principal directions of that deviation matrix.

**P5d (POST-HOC -- see the box above; NOT scored, recorded only so my wrong prior is auditable).**
`||Pi_6 e||^2/||e||^2` is in **[0.15, 0.40]**, i.e. only modestly above the isotropic null
`6/d = 0.154`. Rationale: `pool-error-is-68-percent-common-mode` says most of the error is a shift
every member shares, which is orthogonal to the deviation span by the S23-L9 identity, and nothing
on record suggests the remainder aligns with the pool's own leading modes.
*Falsified* if it exceeds 0.50 -- in which case the subspace is NOT the obstruction and the
combination question is about coefficients, not about span.
