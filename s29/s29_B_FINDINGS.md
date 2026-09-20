# LANE B -- THE COMPATIBILITY HAMILTONIAN, AND THE TAIL-THEN-AGGREGATE OBJECTIVE. SPRINT 29. FINDINGS.

Pre-registration `s29/PREREG_S29_B.md` (base, commit `1c345f07`; addendum 1 `9d745692`;
addendum 2 `d3ccf6b6`) -- every one committed before the number it governs. Code
`s29/s29_B_compat.py` (the three matrices, the two controls, the independent hop-gradient
estimator, the exact ground state, the readouts), `s29/s29_B_tta.py` (the tail-then-aggregate
objective with its envelope gradient, the flatness report, the subset search),
`s29/s29_B_analyse.py` (statistics). Tests `tests/test_s29_B.py` (18 pass). Results
`s29/results/s29_B_*.json|jsonl`. Ledger: S29-L25 (the set-equality counterexample). Tiers as in
S12 to S28: DEMONSTRATED / ORACLE DIAGNOSTIC / DERIVED / REFUTED / OPEN. Every number carries its
artefact path.

BASIS NOTICE. Only the POINT CLOUD appears so far (production 3.048338 at n = 126, 3.252928 on
the 12 trainability targets, `s27/results/s28_B_rows.jsonl :: rmsd_dis75`). No built-chain number
has been produced by this lane; no endpoint has been claimed.

ORACLE NOTICE. Every RMSD below is an ORACLE evaluation of an ACHIEVABLE (native-free) selection.
The matrices, the ground states, the subsets, f, the frame and the readouts never see a native
(`tests/test_s29_B.py :: test_selection_is_native_free_and_nan_poison_is_bit_identical`).

---

## 1. DERIVED -- THE SIGN-MIXING LEMMA, AND WHY A p-READOUT CANNOT CONSUME A CENTERED HAMILTONIAN

Both common-mode-removed compatibility operators annihilate the uniform vector: A_c 1 = H A H 1 = 0
because H 1 = 0, and G 1 = Delta (Delta^T 1)/n_res = 0 because deviations from the mean sum to
zero. Hence every eigenvector with lambda != 0 is orthogonal to 1 and **has entries of both
signs**; the extremal eigenvector is a signed contrast between the two poles of a principal shape
mode; and every deployed readout is a function of p = psi^2, so it is blind to that sign and
averages the two poles back toward the pool mean. Asserted in code
(`tests/test_s29_B.py :: test_sign_mixing_lemma_every_nonzero_eigenvector_has_both_signs`).
Reached independently of, and agreeing exactly with, lane T's S29-L11(d) pole symmetry (my prereg
section 3.2 was written before I read that entry; the entry pre-dates my commit on the clock, so
the claim is independence of reasoning, not of clock).

**Consequence, and it is the closure of the brief's question**: the freedom a centered off-diagonal
term adds lives in the PHASES, and the readout is where it is destroyed. Repairing that needs a
signed readout, and lane O's S29-L21 has now priced the whole signed-readout family's ORACLE
ceiling at **exactly 0.0000 A** (best global eta along the pool's first shape mode is zero; the
per-target sign is positive on 52% of targets, a coin flip; leave-fold-out eta is +0.0071 A
WORSE). So the route is closed at its ceiling, not merely unpromising.

## 2. DERIVED -- THE STABLE-RANK BOUND: NO GRAM OF STRUCTURAL DEVIATIONS CAN BE GRADIENT-VISIBLE HERE

By lane T's law Var[dF/dtheta_0] = r_stable(M)/D^2 at unit spectral norm (S29-L11(b)), "within 30x
of the diagonal terms' 3.051e-2 at D = 512" requires r_stable >= 512^2 x 3.051e-2 / 30 = **266**.
rank(G) <= 3 n_res - 3 = 45 EXACTLY on posed windows (Kabsch removes the three translations
identically; the three rotations only to first order, and in practice exactly: measured rank 36 at
n_res = 14, which is 3 n_res - 6). So r_stable(G) <= 45 < 266 **by construction**, and no
re-weighting of a deviation Gram can fix it. Registered in prereg section 3.3 before measuring and
asserted in `tests/test_s29_B.py :: test_G_is_psd_low_rank_and_shares_the_pool_covariance_spectrum`.
Measured r_stable at n = 9 is far below even that cap: 1.036 (A) / 1.591 (A_c) / 1.675 (G), median
over 12 targets (`s29/results/s29_T_grad_rows.jsonl`, reproduced independently here to 9.9e-15).

## 3. DEMONSTRATED -- THE SET-EQUALITY THEOREM FAILS ON REAL POOLS (S29-L25)

Exhaustive over all C(500, 2) = 124,750 pairs of the deployed pool, 12 targets
(`s29/results/s29_B_tta_subset_rows.jsonl`, job `s26/jobs_done/s29B_tta_subset.json`, exit 0,
92 s, 0.34 GB): the f-optimal 2-subset under f = the shipped distogram risk of the subset AVERAGE
is **not the energy-order prefix on 11 of 12 targets**, and **not the two best by the per-state
criterion on 11 of 12 either** (agreement 1/12). The greedy-plus-local-search m = 5 subset is
non-prefix on **12 of 12**, mean objective gap **+0.1522** over the m = 5 prefix, with supports
reaching ranks 324 and 438 of 500. Lane T's three-state, two-dimensional witness (S29-L17 section
4.5, S29-L18), confirmed on this instrument. This is the first formulation in the project's record
in which "which set" is a genuine optimisation variable rather than a read-out of a sort.

## 4. DEMONSTRATED -- AND THE ESCAPE BUYS NO ACCURACY (S29-L25)

The same subsets' coordinate averages lose to the m = 5 energy prefix on 10 of 12 targets (5/5
folds the same sign, effect +0.2515, 0.52x MDE) and to production on 8 of 12 (+0.4278, 0.59x MDE).
**Every RMSD contrast is NOT MEASURED at n = 12** and is reported as direction only. Lane T's F5c
second clause therefore does not fail: the aggregate-optimal subset is not better than production,
this is not the sprint's opening, and nothing goes to 126. Worst cell 2MP9: the objective moves the
support from the prefix to ranks 44 to 438 and the ORACLE RMSD goes 1.188 -> 3.938, the S29-L2
anti-correlation doing exactly what it says.

## 5. DERIVED AND MEASURED -- THE FLATNESS GATE, AND A CORRECTION TO ITS FRAMING

Registered before measuring (prereg addendum 2, B2.3): the f term of the tail-then-aggregate
objective is flat on **exactly** the CVaR term's flat subspace, because for y above the VaR, p_y
does not appear in R_alpha and perturbing it moves neither the tail set nor the quantile, so
dR/dp_y = 0 by the same envelope argument that gives dCVaR/dp_y = 0. Hence the raw 85.5% **does
not fall**, and the quantity a flatness gate must read is the OVERLAP: among the directions along
which the EMITTED STRUCTURE moves, how many does the objective see. Measured on 1A13
(`s29/results/s29_B_tta_flat_rows.jsonl`): flat_cvar = flat_f = flat_info_bearing = **0.8611**
(440 of 511) identically at every lam in {0, 0.1, 0.3, 1, 3}; the TTA readout's own flat set is the
same 440; **overlap = 0.0000**, i.e. all 71 directions the readout consumes are seen by the
objective, against the deployed readout, which moves along **0 of 511** continuous directions.

## 6. OPEN

- Measurement 1's 12-target independent rows (job `s29B_grad12`) and measurement 2's 12-target
  ORACLE ground-state sweep (job `s29B_gs12`) and the 12-target flatness run (`s29B_tta_flat2`)
  were queued behind `s26/launch_cap.json` max_concurrent 8 at the time of writing.
- Measurement 5's ENDPOINT arm is HELD, on the coordinator's instruction, for lane D's in-band
  ranking experiment: f must clear the meter's numbers 3 and 4 before 126 targets are projected
  through it.
- G_res (deviations residualised on the posterior's own predicted deviation) is registered as NOT
  RUN with its reason (prereg section 6).

## 7. WHAT THIS LANE WOULD SAY IF IT STOPPED HERE

The off-diagonal route on the CANDIDATE-INDEX register is closed, and not for S28's reason. S28
closed it because one similarity measure was near rank one; lane T has now shown that centering
removes that degeneracy and makes the gradient decay WORSE, and that no dense kernel or deviation
Gram can be gradient-visible at D = 512 by a rank argument. What remains true of the class is
structural, not spectral: every common-mode-removed operator's eigenvectors are sign-mixed, and
every deployed readout is sign-blind, so the term's only new content cancels at the readout -- and
the signed readout that would consume it has an ORACLE ceiling of zero (S29-L21). The live object
is the other half of this lane: the tail-then-aggregate objective, which is the first formulation
here where the set is a real variable. It escapes the theorem (S29-L25) and, on the same
measurement, buys no accuracy -- which is the sprint's thesis, that the barrier is information and
not the shape of the optimisation, measured rather than asserted.
