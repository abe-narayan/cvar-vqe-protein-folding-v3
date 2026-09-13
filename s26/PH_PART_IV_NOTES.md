# PART IV NOTES -- THE PHYSICS, FOR A READER WHO KNOWS NO PROTEIN PHYSICS (lane PH, S26)

For the report writer (lane E). Plain sentences, every number with its artefact path. Two pages.
Basis is stated at every RMSD: point cloud (`rmsd_avg`, 3.0483 A), built chain (`rmsd_arm`,
3.2148 A, the production result), relaxed chain (`rmsd_full`, 3.2355 A).

## 1. What the two "physics" scores are

The project has two ways of asking "is this a plausible peptide?" that do not use the native.

**Legacy** is an eleven-term score (`core/energy.py`, `DEFAULT_WEIGHTS`, never fitted): steric
4.0, contact 1.0, local hydrogen bond 1.0, long-range hydrogen bond 3.0, helix cooperativity 2.0,
sheet cooperativity 2.0, solvation 0.5, electrostatic 1.0, aromatic 0.8, torsion 0.15,
compactness 0.4. Each term looks at the backbone and CB atoms only and adds a number; the
weights are hand-set constants. It costs 0.3 ms per structure.

**AMBER** is a real molecular force field, ff14SB with the GBn2 implicit-solvent model, evaluated
by OpenMM (`core/amber.py`). Every heavy atom of every side chain is built onto the backbone by a
deterministic builder (`sidechains.py`, a fixed rotamer, no scan), hydrogens are added, and the
energy is the sum of bond, angle, torsion, non-bonded (Lennard-Jones plus Coulomb) and solvation
terms. A single evaluation costs about 9 ms; a restrained minimisation to convergence costs 6 to
12 s (`s16/repair_FINDINGS.md` section 2.4).

## 2. What "H = E o Relax50" means (s20 L6)

When AMBER is used as a Hamiltonian for the quantum selector, the structure it is handed is an
ideal-geometry chain that has never been relaxed, and its energy there is not a finite number on
42% of the lattice register (`s20/LEDGER.md` L6: finite at Relax_1 112/192, at Relax_50 186/192).
The deployable object is therefore "run 50 steps of restrained minimisation, then read the
energy". That operator, not the raw force field, is the Hamiltonian; the relaxation is what makes
the objective defined, and the Spearman between the energy after 50 steps and after 1 step is
0.36 to 0.88 across three targets, so the ordering is still moving when the cap stops it.

## 3. Legacy is a compactness model; AMBER has the opposite sign

Measured on the shipped K=500 pools of all 126 targets (`s25/results/phys_landscape.json`,
`summary`): the 75 candidates Legacy likes best are 0.758 A MORE compact in radius of gyration
than the pool (SE 0.023) and Legacy's rank correlation with Rg is +0.60. The 75 AMBER likes best
are 1.103 A MORE expanded (SE 0.042), rho(AMBER, Rg) = -0.27. They also disagree with each other
about ordering: rho(Legacy, AMBER) = -0.090 (SE 0.018), reproduced on three instruments
(`s20/results/c_q1.json` -0.0886, s24 -0.0829). Legacy's whole-pool correlation with the true
RMSD is +0.307, but a score that only reads Rg gets +0.279, so almost all of Legacy's apparent
skill is "prefer compact". AMBER's is -0.027 whole-pool and +0.047 with Rg removed.

## 4. Why both rank worse than a random subset

Used as the selector over the pool, with the identical coordinate-average readout and a genuine
CVaR-VQE for all seven configurations (`s25/results/phys_suite.json`, `configs_rank`,
`random_null_rank`): Distogram 3.058, AMBER+Distogram 3.132, Legacy+Distogram 3.215, all three
3.253, random 75-subset 3.425, Legacy+AMBER 3.674, Legacy 3.755, AMBER 3.881 A (point cloud).
Legacy alone is +0.330 A worse than picking 75 at random (1.99x its MDE), AMBER +0.455 (2.42x),
5/5 folds each. Permuting a physics channel while keeping its distribution IMPROVES the endpoint
(Legacy +0.327 worse than its own noise, AMBER +0.471; `s25/agentPHYS_FINDINGS.md` section 1.5).
The reason is section 3: the pool's dominant axis is compactness, the distogram already selects
on it, and each energy pushes along that axis in a direction unrelated to which candidate is
right.

## 5. The steric singularity, and where it lives

An unrelaxed ideal-geometry rebuild of a retrieved window puts atoms on top of each other, and the
Lennard-Jones r^-12 wall turns one overlap into an energy of 1e4 to 1e29 kcal/mol (a relaxed
peptide sits at -1170 to -500). Measured: 58.6% of every pool is above 1e4 kcal/mol; the energies
span 15.3 decades; 97.0% of a pool lands inside |z| < 0.1 of a moment z-score; the ten worst
candidates carry 99.66% of the variance (`s25/results/phys_landscape.json`,
`AMB_frac_absz_lt_0p1`, `AMB_decades`, `AMB_top10_var_share`). In the Pauli basis the same fact
reads: 99.6% of raw AMBER's Walsh variance sits on ten configurations out of 4,096
(`s13/results/walsh_xval.json`, `concentration`, amber median top-10 share 0.996), which makes its
spectrum a delta spike, Binomial(m, 1/2) exactly, and says nothing about the physics. In the
torsion-space Hessian the same fact reads: 7.45% of AMBER's modes carry all of its curvature
(participation ratio 0.0745 against Legacy's 0.4221, 30W/0L; anisotropy 18.8 against 5.8,
0W/30L; condition number three orders larger; `s20/results/c_land_report.txt` section 1). And
when AMBER is minimised from a pool member it moves the chain 0.577 rad per coordinate, five
times Legacy's 0.104, because its first steps are clash relief; 72% of the RMSD damage that
minimisation does is the SIZE of that move, not its direction (`s20/results/c_land_null.json`,
toward-member null +0.444 of +0.620).

**New in S26 (`s26/results/ph_reject_census.json`, `singularity`; ledger L23).** On the
9,450 shipped top-75 rebuilds, 40.7 of every 75 have two heavy atoms closer than 2.0 A when
every side chain is built, but only 2.6 of 75 do on the backbone plus CB (S19 measured 2.66). Of
the 5,057 rebuilds above 1e4 kcal/mol, 96.8% have their closest contact on a side-chain atom
(bb-sc 2,627, sc-sc 2,267, bb-bb 163). Within a top-75 the AMBER single point is the minimum
heavy-atom distance to Spearman -0.74. So the singularity is mostly the builder's: the fixed
rotamer that places each side chain onto a backbone it did not see. The backbones themselves are
almost always physically possible.

## 6. The non-monotone standardisation trap

Turning an energy into a score by subtracting its mean and dividing by its standard deviation
looks harmless, but with one 1e28 outlier setting the standard deviation, every candidate below
about 1e12 maps to the same double-precision number. On 40 of 126 targets the moment z-score
CHANGES AMBER's own ordering, with exact tie blocks up to 462 of 500; `np.argsort` then returns
those tied candidates in array order, which is the BLOSUM retrieval order, and that order is not
neutral (rho with the true RMSD +0.054, seven standard errors from zero). Every Angstrom of the
moment-z "advantage" on the AMBER configuration lived on exactly those 40 targets
(`s25/agentPHYS_FINDINGS.md` section 3; `s25/results/phys_suite.json`, `normalisation_fork`).
Rank standardisation is monotone and has no such failure; it is what production uses
(`core.pipeline._zrank`).

## 7. The relaxation, on the production chain (native-free part, S26)

The last production stage relaxes the built chain with a restraint of 10 kcal/mol/A^2 on N, CA
and C. Measured on the 126 emissions (`s26/results/ph_c3_nativefree.json`; ledger L24): the
built chain's own AMBER energy is above 1e4 kcal/mol on 58.7% of targets before relaxation
(median 8.6e4), the relaxation brings it to -560 on average and converges on 125 of 126 (9KAR
ends at +1262), moves the CA trace 0.220 A RMS, and stretches the virtual CA-CA bond from 3.804
to 3.867 A on average, breaking it beyond 4.0 A on two targets (2BP4 5.38 A, 9KAR 4.86 A). What
that step does to accuracy against a random move of the same size is C3 (`s26/C3_RESULT.md`,
after the gate); on the record so far the step costs +0.0207 A [+0.0143, +0.0276] and S16 found a
random displacement of matched size at least as accurate (`s16/LEDGER.md` L27).

## 8. The sentence for the presenter

The two physics scores are real, they disagree with each other about what a good peptide looks
like (one likes compact, one likes spread out), and on this pool both are worse than choosing at
random because neither is looking at the thing that decides accuracy. The all-atom force field's
number on an unrelaxed candidate is almost entirely the distance between its two closest atoms,
and that distance is set by how the side chains were placed, not by the backbone. It belongs in
the pipeline as a validity check, which is what it is, and not as a judge of which candidate is
right, which it is not.
