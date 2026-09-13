# IDEA -- The cis-peptide gap: census, floor, and a two-bond-length projection design (lane PH; mandatory direction)

## Hypothesis
The production projection can only emit trans peptide bonds (`core.geometry.OMEGA_TRANS = pi`;
`core/project.py` `frames` takes a scalar omega; constant 3.804 A virtual bond,
`s26/results/e_reproduce.json`). A native with a cis bond (CA-CA ~2.9 A) therefore has a
representation floor the projection cannot go below. The hypothesis has two parts: (1) the
floor is material on the dev instrument; (2) a per-bond omega (two values, trans and cis) chosen
by a native-free rule (proline-preceding positions) would recover it.

## Why the record does not already close it
The state brief lists the constant virtual bond as a declared defect (section 3) but no sprint
measured its cost. `core/data.py:406-407` drops peptides with any consecutive CA-CA step below
3.5 A from the database, and `core/data.py:697-698` drops fragment windows with such a step, so
the cost on the instrument may be zero BY CONSTRUCTION of the database rather than by the
projection; that has never been stated with a number either. S23 L11 measured cis fractions of
0.51 to 0.64 appearing when only CA is restrained in AMBER (the relaxation manufactures cis
bonds), and S19/S21 used `cis_frac` as a validity axis (+0.315 to +0.424 disqualified two arms,
`s21/LEDGER.md` L36), so cis bonds are on the record only as a repair artefact, never as a
target property.

## Exact falsifier
Part 1 (census, allowed before the gate): the number of dev natives (model 1) with |omega| < 30
degrees, cross-checked by CA-CA < 3.3 A, with the residue after the bond, plus every deposited
model and the pool windows. Registered prediction: zero on model 1 and zero in the pool because
of the step gates; some in other ensemble models; a non-planarity tail. Part 2 (after the gate):
the CA floor of ideal trans geometry on the native's own torsions (`core.data.Peptide.rebuild`
on N/CA/C, verified from source, and the CA-only twin), on cis targets against non-cis, against
the production built-chain cost `rmsd_arm - rmsd_avg`. The design is worth writing iff at least
one cis target exists and its CA floor exceeds the MDE of the chain-cost comparison. Full
pre-registration: `s26/PREREG_cis.md`.

## Expected effect against the MDE
On the 126: zero cis targets, so the two-bond-length projection is worth 0.000 A on the
instrument and the design is written for the world supply (16 containment-fresh targets, 10
amyloid). The non-planarity floor is expected small against the 0.166 A chain cost (S16 measured
the projection's cost as a displacement-magnitude effect, not a representation effect). If the
floor is material without cis bonds, the mechanism is omega non-planarity and the design
generalises to a continuous per-bond omega.

## Cost
Part 1: 2 min CPU, < 0.5 GB. Part 2: 1 min. Design: 1 agent-hour of writing, no implementation
on the production path. Agent-hours total: 2.

## Information value
Turns a declared defect into a priced one; either "zero on this instrument by construction, X on
the world supply" or a real number. Low plausibility of a gain, high value as a closed line in
the report.
