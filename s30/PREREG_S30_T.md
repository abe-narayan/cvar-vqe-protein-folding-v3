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
