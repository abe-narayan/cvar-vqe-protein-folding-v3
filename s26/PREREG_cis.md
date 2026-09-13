# PREREG -- THE CIS-PEPTIDE GAP (lane PH, S26)

Written 2026-09-13 00:40, before `s26/ph_cis.py` produced a number. Never edited after; addenda
appended. Mandatory tournament direction. Code: `s26/ph_cis.py`. Results:
`s26/results/ph_cis_census.json` (part 1), `s26/results/ph_cis_floor.json` (part 2).

## 1. The gap

The production projection builds every chain at ideal TRANS geometry: `core.geometry.OMEGA_TRANS
= pi`, a scalar omega in `core.project.frames` (line 305) and `build_ca` (line 428), a constant
virtual CA-CA bond of 3.804 A (`s26/results/e_reproduce.json`, `chain_mean_bond` 3.80395,
`chain_bond_sd_max` 1.8e-15). A cis peptide bond has CA-CA about 2.9 A and cannot be represented.
The state brief lists this as a declared defect (section 3). What is not on the record is how
much it costs on the instrument, and whether it can cost anything at all given the database
filters.

## 2. Part 1, the census (allowed before the gate, ledger L5)

For each of the 126 dev targets: parse N, CA, C of model 1 (the model the instrument scores
against, `core.geometry.parse_pdb(model_index=0)`), compute omega_i = dihedral(CA_i, C_i,
N_{i+1}, CA_{i+1}); cis iff |omega| < 30 degrees; cross-check with consecutive CA-CA < 3.3 A;
record the residue after the bond (proline or not); record the omega deviation |180 - |omega||
per bond. Separately, the same over every deposited model. For the pool: consecutive CA-CA <
3.3 A over every window of the universe (`W`, CA only), the K=500 pool and the production
top-75. Never enumerates any directory; iterates over the 126 dev ids only; reads no `rr`, no
`nat_ca`, computes no RMSD.

**Registered expectation, and why it is a prediction about the FILTER rather than about
peptides.** `core/data.py:406-407` drops any peptide whose consecutive CA-CA step is below 3.5 A
or above 4.1 A, and `core/data.py:697-698` drops any fragment window containing such a step. The
126 dev targets are drawn from that database (`s7/debias.py:103-120` iterates `db.load()`), and
every pool window passes the fragment gate. I therefore predict: zero model-1 natives with a cis
bond by either criterion, zero cis windows in the universe, and a universe minimum step above
3.5 A. I predict non-zero cis bonds in OTHER ensemble models on some targets (NMR ensembles are
not filtered model by model), and a non-planarity tail (omega deviation > 20 degrees) on a
minority of bonds, which is what the projection's fixed omega actually fails to represent on
this instrument.

Falsifier of that reading: any model-1 cis bond found. Then the filter is not doing what the
code says, and the floor on that target is the part-2 measurement of record.

## 3. Part 2, the projection floor (after the gate)

`core.data.Peptide.rebuild` was read from `core/data.py:410-419` before use: it is the RMSD
between the native N/CA/C backbone and the chain rebuilt from the native's own phi/psi through
`core.geometry.build_backbone` (ideal trans geometry) after Kabsch superposition. That is the
representation floor of the production projection on the native's own torsions, on the N/CA/C
basis, and every database entry has it below `REBUILD_TOL = 1.5`. Part 2 also computes the
CA-only floor (the endpoint's basis), and reads the production `rmsd_arm - rmsd_avg` per target
(the built-chain cost, where a wrong virtual bond would show).

ORACLE DIAGNOSTIC throughout (reads the native). Compared: cis targets vs non-cis targets (if
any cis target exists); the CA floor against the chain cost, paired, `ST.compare`, fold CI, MDE;
Spearman of the CA floor and of the chain cost against the per-target maximum omega deviation.

Falsifier for "the two-bond-length projection is worth designing on this instrument": the CA
floor on cis targets exceeds the MDE of the chain-cost comparison and at least one cis target
exists. If no cis target exists, the cost on the instrument is zero by construction of the
database, the design is written for the world supply only (the 16 containment-fresh targets, 10
amyloid), and the non-planarity floor is reported as the residual representation cost.

Expected: the CA floor is a small fraction of the 0.166 A chain cost (which S16 measured as a
magnitude effect of the projection's move, not a representation effect), and no cis target
exists. If the floor turns out material without cis bonds, omega non-planarity is the mechanism
and the design generalises to a per-bond omega rather than to a two-valued one.

## 4. The design deliverable (not an implementation)

If the cost is material, `s26/agentPH_FINDINGS.md` carries: where omega enters
(`core/project.py` `frames` takes a scalar omega and fills `t[:, 1::3] = omega`; `build_ca` and
`build_ca_exact` pass it through; `core.geometry.build_backbone(phi, psi, omega=OMEGA_TRANS)`),
what a per-bond omega vector changes (the transform at step 3i+1 for each bond; the analytic
gradient is unchanged because omega is not optimised), what native-free signal would choose cis
(only the proline-preceding position is defensible; base rate to be stated from the census), and
the cost of a test. Nothing on the production path is edited.

## 5. Cost and memory

Part 1: ~2 min, < 0.5 GB (parses 126 PDB files and their ensembles through the cached column
reader; loads each `W` bank once). Part 2: ~1 min. Tag CPU. Agent-hours: 1.5.

---
## ADDENDUM 1 (2026-09-13 08:41) -- THE CENSUS RESULT AND THE RE-SCOPE OF PART 2

`s26/results/ph_cis_census.json`, ledger L22. The registered prediction held on model 1 and on
the pool (0/126, 0/2,352,893, universe minimum step 3.5045 A, which is the gate) and FAILED on
the ensembles: 0 of 1,966 deposited models carry a cis bond, where I predicted some. Written
before part 2 runs: the cis-vs-non-cis comparison is empty; part 2 measures the ideal-trans
representation floor (CA and N/CA/C bases) on all 126 as the residual cost of the constant
omega, its Spearman with the per-target maximum omega deviation, and reports 9UV5 (the one
target with a bond beyond 30 deg) by name, never as a subgroup. The falsifier for "the design is
worth writing on this instrument" cannot fire, so the design is written for the world supply
with the base rate stated as 0 of 1,507 bonds here and a literature value with its source.
