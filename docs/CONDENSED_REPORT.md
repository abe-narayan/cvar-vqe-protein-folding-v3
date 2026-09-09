# Where the accuracy isn't - condensed record, Sprints 5-13

Peptides of 9-16 residues. Metric: mean CA-RMSD to deposited model 1.
Instrument: 126 cluster-disjoint tuning targets. The 60-target benchmark was read
once, on pre-registered constants, and has not been touched since.

Web version: https://claude.ai/code/artifact/53485906-4307-4185-b2ff-47dd7673d509

---

## The bottom line

**There is no validated accuracy improvement over the shipped baseline.** On the
held-out benchmark the synthesis architecture scored 2.9610 against the baseline's
2.9507 - paired +0.0103, 95% CI [-0.1596, +0.1803], 31W/29L. The mean transferred;
the effect did not.

Two candidate architectures have been measured to their ceilings.

- The **retrieval pipeline** is signal-limited at selection: a candidate set holding
  a 1.71 A answer is reliably converted into a 3.21 A one, and the last untested
  class of ranker returns a flat learning curve.
- The **torsion pipeline** has a better representation and no usable objective: a
  space holding a 1.59 A answer in ~24 qubits, in which neither energy ranks the
  native above the 32nd percentile.

What survives is not the folding. It is a measurement of how a molecular energy
model's structure determines the trainability of a variational quantum algorithm,
closed end to end with no free parameter, on a genuine force field.

---

## The ladder (126 targets, same metric)

ORACLE arms read the native and are diagnostics, not results.

| arm | mean CA-RMSD |
|---|---|
| annealing on the Legacy energy | 4.624 |
| random sampling, 5,000 draws | 4.522 |
| a constant alpha-helix, zero information | 4.065 |
| annealing on a 1-local torsion prior | 3.969 |
| sequence-only torsion predictor | 3.770 |
| **shipped retrieval pipeline** | **3.213** |
| ORACLE perfect distance objective, m=75 | 2.395 |
| ORACLE torsion space, random-start search | 1.982 |
| ORACLE best member of the retrieved pool | 1.711 |
| ORACLE torsion space ceiling, 24 qubits | 1.594 |
| ORACLE torsion restraints at sigma = 12 deg | 1.486 |

---

## Line one: retrieval, and why selection closed

**The last untested ranker is built and the learning curve is flat.** A
set-transformer over the signed 75 x n_pairs x 18 deviation tensor, 18 hardest
targets held out: emits 3.184 vs production 3.203, d = -0.019 [-0.058, +0.020],
63W/63L. Learning curve n=8 -> n=75 is 3.043 -> 3.026. The identical harness with
one leaked RMSD label reaches 2.534 at n=8. **Signal-limited, not sample-limited.**

**The terminal operator consumes the set mean, not the set best.**
`d_out = 1.16*d_set_mean + 0.04*d_set_best`, R^2 = 0.893, over 882 perturbations.
A perfect rank-1 decision is worth -1.743 A through argmin and -0.029 A through the
m=75 average the system ships. Optimal set size shrinks 500 -> 75 -> 20 -> 5 as the
objective improves.

**Sequence conditioning is worth 0.776 A overall and is NEGATIVE where it is needed.**
Blind pipeline vs shipped: 3.749 vs 2.745 on the 108 ordinary targets (-1.004);
5.425 vs 6.019 on the 18 hardest (+0.594). A pipeline told nothing about the target
sequence beats the real one on the failure class. Mechanism: identity predicts
structure at rho ~ -0.25 among isolated peptides but only -0.03 among protein
fragments, which supply ~80% of every pool. Three independently built routers were
null.

**Ten of the eighteen failures are outside the representation.** Fibril segments and
lasso peptides, 10/18 vs 6/108, Fisher p = 1.2e-6, surviving Bonferroni over every
flag searched and leave-one-out. All six lasso peptides emit >= 6.10 A; three have
zero windows within 2 A anywhere in a 10,000-window universe.

---

## Line two: torsions, and why the objective closed

**The representation is excellent; the library is mostly not the reason.** At k=4
(~24 live qubits) the space contains a 1.594 A answer. But a wrong target's library
costs only +0.135 A, a class-blind single codebook costs +0.135 A, and
uniform-random torsion states still reach 2.698 A. Of the 1.104 A the library buys
over an information-free space, **88% is generic Ramachandran and 12% is anything
sequence-related**; 0.388 A of the ceiling is the privileged oracle start.

**Neither energy ranks the native.** Legacy's rank correlation inside the low-energy
decile is +0.043; raw AMBER -0.088. The native sits at the 32nd-40th percentile of
both. A monotone log compression leaves rho unchanged at -0.036 while collapsing the
range from 16.1 decades to 1.8 - **the damage is in the ordering, not the scale**.
On nine fully enumerated targets (262,144 configs each), Legacy's certified global
optimum is **+0.139 A worse than random sampling**.

**Optimising harder makes the structure worse.** 3.764 at 10 evaluations -> 3.667 at
300 -> 3.920 at the certified global optimum. The turn is governed by the fraction of
space explored, not the evaluation count. A VQE at n=12-16 with 10^4 evaluations sits
on the improving limb of a curve whose limit is known to be bad. **It will look like
it is working.**

**The sequence carries no phi.** Leave-fold-out mean absolute torsion error, degrees:

| model | phi | psi |
|---|---|---|
| sequence-blind corpus marginal | 36.4 | 72.8 |
| residue-class prior, no learning | 36.9 | 68.9 |
| full sequence context and properties | 36.1 | 62.4 |

A model seeing the entire 15-residue context predicts phi no better than one seeing
no sequence at all. The whole measurable channel is 10.4 deg of psi. Best
sequence-only builder: sigma 67.7 deg, 3.770 A emitted, and a perfect confidence gate
still loses at 3.478. The incumbent pipeline corresponds to sigma ~29 deg, which is
the number to beat merely to tie.

---

## What survives: the trainability result

**An exact locality theorem in torsion space.** Under an ideal-geometry backbone
builder, `d_ij` depends on exactly the `j-i-1` residues strictly between i and j,
contiguous; agreement 1.0000, zero counterexamples at machine precision, holding at
termini, on GLY/PRO, and under cis-omega. All-atom follows four analogous exact rules
(CA `i<m<j`, N `i<=m<j`, C/O `i<m<=j`, CB `i<=m<=j`).

Consequences: a separation-8 pair is a 14-qubit interaction at k=4, so **no 2-local
Ising form of a distance-based molecular objective exists in this encoding**; and
both energies are full-register, so "AMBER is less local than Legacy" is a category
error - recorded as refuted.

**The spectrum of a force field predicts its gradient variance, with no free
parameter.**

| energy model | mean Pauli weight | measured / predicted |
|---|---|---|
| Legacy, coarse-grained | 2.236 | 1.006 |
| AMBER ff14SB/GBn2, conditioned | 3.015 | 1.001 |

AMBER is higher on 79/79 cells. The only approximation is dropped cross-covariances
between Pauli strings, and the ratio is that approximation's test. The first attempt
was an artefact: raw AMBER's top-10 configurations of 4,096 carried a median 99.6%
of its Walsh variance, and a constant-plus-spike has weight spectrum exactly
Binomial(m, 1/2) - measured 6.001 against predicted 6.001. Ninety-ninth-percentile
winsorisation was insufficient; only monotone rank-preserving conditioning works.
Per physical term, AMBER's non-bonded contribution has covariance share 1.000 and
Legacy's steric term 0.955: both models are a steric potential plus rounding error.

**Three refutations worth keeping.**

1. **Cost-locality does not explain trainability in this regime.** The measured
   ansatz kernel is flat in Pauli weight; a maximally global observable's gradient
   variance is within 1.2-4.8x of a single-qubit one. What matters is n. The standard
   theorem needs local 2-design blocks this ansatz does not have, so it licenses no
   prediction at 6-18 qubits.
2. **QNG has nothing to fix.** The metric is full rank at every parameter, depth and
   size, with `g_ii = 0.2500` exactly, off-diagonal correlations shrinking with n,
   and exactly `I/4` at depth 1.
3. **The metric contains no Hamiltonian.** Bit-identical across energy models at
   matched parameters (0.000e+00). The energy model does not reshape the manifold; it
   selects which region the optimiser visits. Now a passing unit test.

**Two more, on discipline.** A zero-information constant alpha-helix (phi=-63,
psi=-42) beats the mandatory random control by 0.457 A [-0.822, -0.090], so "beats
random" is not evidence of anything in this representation. And SPSA optimises the
AMBER objective best of five arms while returning the worst structure (+0.333 A vs
random): objective quality and structural quality are separate axes.

---

## The corrections ledger (the ten that changed a decision)

| was | now |
|---|---|
| The pool caps accuracy at 2.406 A | Measured through a weak terminal; averaging the same pool gives 1.925 |
| No objective upgrade reaches 2.0 A | True only at fixed cardinality; with a co-optimised terminal, ~1.99 A |
| AMBER is a less local observable than Legacy | Category error; both full-register, difference smaller than reported |
| AMBER's weight exceeds Legacy's on 14/14 cells | Measuring one steric clash; a delta spike is global by arithmetic |
| Optimising Legacy is worse than not optimising | Reverses under stronger optimisers; say energy and structure are uncorrelated (rho -0.006) |
| The torsion library's sequence conditioning is the asset | 88% of it is generic Ramachandran |
| Clustered restraint gaps are the worst case | Opposite; terminal dropout is 0.40-0.50 A cheaper than uniform |
| A sequence-only torsion route lands at 2.4-2.9 A | Measured 3.770 A; the published read was MSA-mode, not single-sequence |
| A 0.6 containment threshold screens leakage | At the null (0.56-0.63 for random matched sequences); 4 targets contain themselves |
| An AMBER single-point costs 6 ms | 28 ms; the first timing re-evaluated the same state and hit a result memo |

Twenty-three corrections stand across the full record; these are the ones that moved
a decision.

---

## Where this leaves the work

**Closed - do not re-fund.** Further in-band rankers over retrieval pools; post-hoc
correction of the distance objective in any form (two independent deaths, different
mechanisms: calibration, then error coherence); growing the fragment library (~0.43 A
per decade, 800x for 0.8 A); multi-piece assembly as a generation lever,
ideal-template banks, finer structural alphabets, key fusion, more protein-fragment
training data; putting either energy in a search objective, whose honest role is a
validity filter.

**Open, in order of value.**

1. **A many-body objective trained against accuracy.** The enumerated ground truth
   makes this supervised, not physics: 2.36 M labelled structures with true RMSD,
   eleven energy terms and a prior. The prize is the 0.95 A between the 1-local
   optimum (1.921) and the space best (0.969).
2. **Chemical-shift-derived torsion restraints**, coverage measured from BMRB
   deposits first, before any structure code is written. Reaches 1.486 A at full
   coverage and sigma 12 deg; needs >= 90% coverage; the coverage model was just
   corrected in its favour. Carry the caveat: deposited coordinates were solved using
   those shifts, so any such arm is NMR-restrained prediction in its own column.
3. **Publishing the trainability half.** Complete, corrected, cross-validated, and
   novel in the one dimension the literature leaves open - varying the energy model
   as the independent variable under matched variational conditions.

**The validation constraint.** No adequately powered fresh benchmark can be built.
The three existing instruments spend all 204 identity clusters of 9-16-residue
peptides in the corpus. Outside it: 700 single-chain entries of that length, 28 from
this year, 146 absent from the corpus entirely. Built out in full, those 146 yield
**16 usable targets after quality and leakage gates, 10 of them amyloid fibril
segments** - the class already excluded as not folded in isolation. Any future claim
must name which of three imperfect instruments it rests on.

---

Every number here is a paired comparison on the 126-target tuning instrument with a
bootstrap interval, a win/loss split, per-fold values and a drop-top-ten check,
unless marked ORACLE. No architecture decision was made on the benchmark.

Full record: `FINDINGS.md`, `s12/SPRINT12_DOSSIER.md`, `s13/SPRINT13_DOSSIER.md`,
`s13/coord_FINDINGS.md`, twenty findings files, seven figures in `s13/figures/`.
