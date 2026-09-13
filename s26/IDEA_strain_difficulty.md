# IDEA -- The relaxation's own strain as a native-free difficulty signal (lane PH; own idea 2)

## Hypothesis
The production relaxation reports, for every target and at no extra cost, the built chain's own
AMBER energy before relaxation (`amber_e0`, the strain of the emitted chain: 6.0e5 kcal/mol on
1A13), the energy after (`amber_e1`), how far the restrained atoms had to move (`amber_moved`)
and the residual bond+angle strain (`amber_strain_after`), all in
`bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json`. The hypothesis is that these native-free
scalars predict the per-target error of the built chain (`rmsd_arm`): a chain that packs itself
into a hard clash, or that the force field has to move far to make physical, is a chain the
averaging produced from an inconsistent set. This is NOT an accuracy lever; it is a calibration
output ("the physics reports when the answer is untrustworthy") that the presentation currently
lacks, and it would be the first use of the physics that the record does not already price at
zero.

## Why the record does not already close it
- S10 ("recognition is the barrier") and the per-target routers of S22 L7 / S23 L7 fitted
  native-free features to predict the optimal set size or scale and failed held-out. They never
  used the relaxation's energies, which exist only after the emission.
- S20 L4: circuit-side landscape metrics are difficulty proxies. That was the VQE's landscape,
  not the emitted chain's strain.
- `in-band-ordering-is-per-target` (memory): native-free compactness proxies reach 0.24 to 0.37
  against the oracle's 0.909 for the per-target SIGN; strain was not among them.
- S20 L12: AMBER's `spec_skew` at the start predicts the RMSD its minimiser reaches (rho -0.487,
  n=30). That is the same mechanism seen in torsion space on 30 targets; it has never been
  measured on the production emission at n=126.

## Exact falsifier
Spearman across the 126 targets between each of {log e0, log(e0 - e1), moved, strain_after} and
`rmsd_arm` (ORACLE DIAGNOSTIC evaluation of native-free quantities), with a fold-clustered
bootstrap CI; a label-permutation null. The signal exists iff at least one scalar reaches |rho|
>= 0.25 with the CI excluding zero and the sign replicating on 5/5 folds. Also the practical
form: does the top quartile by strain contain more of the FAIL18 targets than chance (Fisher
exact, one-sided, alpha 0.05, Bonferroni over the four scalars)? Confound control: chain length
n and Rg of the chain are partialled out (longer chains have more atoms and more strain), and
the result is reported with and without the partial.

## Expected effect against the MDE
Unknown; a correlation of 0.2 to 0.3 is my guess, which is the size at which the routers failed
to convert a signal into an operator. The difference here is that nothing is converted: the
deliverable is a confidence flag beside the emitted structure, judged by its calibration curve,
not by an RMSD gain. Plausibility 0.5 for the signal, 0 for an RMSD gain.

## Cost
Zero compute beyond reading 126 cached records: 1 min CPU, < 0.3 GB. Agent-hours: 1.5.

## Information value
Moderate. A positive gives the presentation a defensible use of the physics ("validity plus a
warning"); a negative closes the cheapest remaining difficulty predictor. It also feeds C3: if
the strain does not predict the error, the relaxation cannot even tell which targets it is
about to make worse.
