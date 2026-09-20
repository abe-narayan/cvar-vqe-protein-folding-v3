# L_5 -- PEPTIDE STRUCTURE FROM SEQUENCE AT 9 TO 16 RESIDUES: THE PUBLISHED CEILING (brief topic 5)

Lane L, 2026-09-20. The question: what accuracy does the field reach at this length, which of
those inputs is legitimately available here (CPU only, no structure-trained predictor because of
leakage against the natives), and how does the project's 3.21 A compare?

Answer, stated first. The published native-free numbers at 9 to 25 residues cluster at 2.0 to
2.6 A, and the project's 3.2126 A built chain is worse than all of them -- but NOT ONE of those
numbers is measured on a comparable instrument, and three specific differences account for the
gap before any question of method quality arises: (1) benchmark composition (25 to 56 curated NMR
peptides, mostly with regular secondary structure, versus 126 identity-clustered PDB targets that
include the amyloid and lasso classes the project's FAIL18 is made of); (2) best-of-N versus
single-answer reporting (several published figures are "near-native in the top 5", and AF2's own
rank-1 model is 0.2 to 1.1 A worse than its best-of-5); (3) every one of those methods generates
with the same energy or restraints it selects with, which the project cannot do with a retrieval
pool. A like-for-like comparison does not exist in the literature, and manufacturing one is not
possible without re-running those methods on this benchmark, which the leakage rules and the
hardware both forbid.

---

## 1. The published numbers, with what each one actually measures

| method | benchmark | length | reported | what the number is |
|---|---|---|---|---|
| PEP-FOLD 1 (Maupetit, Derreumaux, Tuffery, NAR 37:W498, 2009) | 25 NMR peptides | 9-25 aa | "the PEP-FOLD LEC reproduces the NMR structure at 2.6 A cRMSD", averaged over 25 peptides and 5 runs | the LOWEST-sOPEP-ENERGY cluster centroid: a genuine native-free single answer |
| PEP-FOLD 2 (Shen et al., JCTC 10:4745, 2014) | 56 peptides | 25-52 aa | near-native for 95% of targets (Rosetta 88%) | a hit rate over the generated set, not a selected-answer RMSD; and out of our length range |
| PEP-FOLD 3 (Lamiable et al., NAR 44:W449, 2016) | 56 peptides | 25-52 aa | "a near-native or native conformation in the top five best scored models for 80% of the targets" | BEST-OF-5 after sOPEP/Apollo ranking |
| PEP-FOLD 4 (Rey et al., NAR 51:W432, 2023) | 588 peptides (10-40), plus small pH sets | 10-40 aa | no explicit RMSD table; "the three methods performed similarly on a total 17 peptides" vs trRosetta and AF2 | qualitative; its advantage is pH-dependent and poly-charged cases |
| APPTEST (Timmons & Hewage, Brief Bioinform 22:bbab308, 2021) | 42 short peptides | 9-25 aa | 1.96 A (PEP-FOLD 2.05, PEPstrMOD 4.66) | selection by lowest XPLOR-NIH energy / CYANA target function -- native-free and deployable |
| AlphaFold2 on peptides (McDonald et al., Structure 31:111, 2023) | 588 NMR peptides | 10-40 aa | best-of-5 by class: disulfide-rich 2.2, helical-membrane 2.3, beta-hairpin 2.9, mixed-soluble 4.4, helical-soluble 4.5 A (preprint, 155-peptide version) | BEST-OF-5 against the NMR ensemble; rank-1 is 0.2-1.1 A worse; "no correlation between the first three ranks assigned by AF2 and the structure that gave the lowest Ca RMSD" |
| MD (Lindorff-Larsen et al., Science 334:517, 2011) | 12 fast-folding proteins (chignolin 10 aa, Trp-cage 20 aa) | 10-80 aa | folds to the experimental native repeatedly | an equilibrium population at 100 us to 1 ms per system |

This project: 3.2126 A mean built-chain CA-RMSD, 126 targets of 9-16 residues, single deployable
answer, no best-of-N (`s27/results/chain_rows.jsonl :: DIS`; point cloud 3.0483).

## 2. Why the comparison is not like-for-like, in three specific ways

### 2.1 Benchmark composition

The published sets are curated NMR peptides in aqueous solution with regular secondary structure.
The project's instrument is 126 identity-clustered targets drawn from the PDB, and its hard
stratum is chemically identifiable: S12 measured that "56% of FAIL18 is two chemically
identifiable classes (steric-zipper amyloid segments, lasso peptides) that no linear-window
retrieval can ever represent" (`s12/lit_FINDINGS.md`). PEP-FOLD states its own scope excludes
peptides "not bound to membrane or ligands, and not complexed to or stabilized by metal ions".
These sets exclude, by construction, much of what makes this benchmark hard. Any comparison that
ignores this is comparing difficulty, not method.

### 2.2 Best-of-N versus a single answer

PEP-FOLD3's 80% and McDonald's per-class numbers are best-of-5 quantities. The project's record
is unusually clear-eyed about what that costs: the same paper that supplies the AF2 numbers also
shows the rank-1 model is 0.2 to 1.1 A worse than best-of-5, with no rank-RMSD correlation at
all (L_1 section 4.1). The project reports a single deployable answer and would have to be given
the same best-of-N allowance to be compared; its ORACLE top-75 ceiling of 2.31 A and pool-best of
1.71 A (charter finding 10) are the corresponding best-of-N quantities, and they sit inside the
published band. That is the fair reading: the project's GENERATION is competitive with the
published field at this length; its SELECTION is what is missing, which is finding 8 again.

### 2.3 Generation and selection share an energy in every published method

PEP-FOLD assembles fragments under sOPEP and then ranks clusters by sOPEP; APPTEST folds under
predicted restraints and then ranks by restraint violation; AF2's rank is its own pLDDT head. In
each case the selector scores an ensemble it (or its own model) produced, so the score is a
convergence diagnostic over a self-consistent set. The project retrieves 500 real windows from
other proteins and scores them with an independently constructed objective. L_1's conclusion
applies directly: the only native-free selectors with peptide-length evidence work in the
self-consistent regime, and the one that does not (AF2's pLDDT ranking its own five models) has
no in-band skill.

## 3. Which of these inputs is legitimately available here

| input | available? | why |
|---|---|---|
| a 27-state structural-alphabet profile predicted from a PSSM (PEP-FOLD's retrieval key) | NO -- CLOSED, and I nearly re-imported it | See section 3.1. |
| a torsion head (phi/psi predicted directly; APPTEST) | NO | S13: phi carries no sequence signal at peptide length (36.1 deg with full context vs 36.4 sequence-blind); direct build emits 4.151 A |
| a fold-from-restraints terminal operator (APPTEST, AlphaFold1) | YES in principle | `core/project.py` already does L-BFGS over (phi,psi) with multi-start; what is missing is running it against the distogram rather than against a retrieved coordinate average (S12 Part 4). An open architecture question for lanes M/X, not a QA import. |
| a structure-trained predictor (AF2, ESMFold, OmegaFold, trRosetta) | NO | leakage against the natives, and ESMFold is infeasible on this box (S26 VII.2 / L13) |
| sOPEP or another peptide-tuned coarse-grained force field | NO VALUE | same functional class as Legacy's 11 terms, which is measurably worse than a random subset as a ranker on this pool (S25 L16) |
| converged MD free energy | NO | 10^5-10^6 CPU-hours per peptide against 6.43 core-equivalents |
| chemical-shift torsion restraints (CS-Rosetta/TALOS-N route) | CLOSED | 54/126 coverage and ORACLE-perfect torsions still leave the instrument at 2.021 A (memory: torsion-restraints-reach-the-target) |

### 3.1 The retrieval key: closed, and a near-miss by this lane

PEP-FOLD's architectural claim is that the right fragment key is PREDICTED LOCAL CONFORMATION,
not sequence similarity, and S12 recorded a striking pool-best diagnostic for it (2.284 -> 1.640 A
on FAIL18, near-native recall 16/500 -> 92/500, robust to 30% prediction error). I initially wrote
this up as the strongest unexploited lead in the peptide literature. THAT IS WRONG, and the record
says so in two places I had to go and check:

- S13 built the predictor: "a leave-fold-out 4-state torsion-bin predictor, 0.690 accuracy
  overall, 0.517 on FAIL18 (majority baseline 0.562) -- it failed as a retrieval key"
  (`s13/BRIEF.md`). The ABEGO-style key was tried and did not work.
- The 22-key screen over the whole window universe plus 24 chain arms on 126 paired targets
  (`docs/FINDINGS.md`) tested every structural key available: `disto` (the shipped distogram as a
  key), `dconf`, `dshort` (a secondary-structure key, which is the closest thing to a structural
  alphabet), and fusions. The distogram key produces the largest pool-DISTRIBUTION move in the
  project (best-window percentile 39.3 -> 16.5 median, sub-2 A count 54.5 -> 101.5, pool mean
  4.453 -> 3.605) and yet makes the pool WORSE where it matters: pool-best 1.711 -> 2.161 and the
  distance-ORACLE ceiling 1.994 -> 2.542, because "the distogram key concentrates the pool on one
  predicted structural type: many near-native candidates on the targets it gets right, none at
  all on the targets it does not". On the chain, no key beat BLOSUM (incumbent 3.454; best
  alternative `fuse_be` at 3.447, CI [-0.100, +0.083]). S26's report lists the retrieval key among
  the closed items.

So the published architecture's central idea has been tested here, twice, in its own terms, and
it fails for a stated mechanism (concentration on the predicted type, which is the typicality
axis again). Recorded as a lane-L near-miss because it is exactly the failure mode the brief
warns about -- importing something back that the record closed -- and the S12 line that tempted me
is a FAIL18 pool-best diagnostic, not a chain result.

## 4. So how does 3.21 A compare?

Stated honestly, three ways, because one number would mislead:

1. Against published SELECTED answers at 9-25 aa (PEP-FOLD 2.6, APPTEST 1.96): the project is
   0.6 to 1.25 A worse, on a harder and larger benchmark, with a selector that does not share an
   energy with its generator. The comparison is unfavourable and should be reported as such,
   with the three caveats of section 2 attached, not used to explain it away.
2. Against published BEST-OF-N at 10-40 aa (AF2 2.2-4.5 by class): the project's matched
   quantities are ORACLE top-75 2.31 A and pool-best 1.71 A, which are inside or below that
   band. Generation at this length is not the project's problem.
3. Against the only zero-information reference that exists here: production 3.2126 vs the
   random-75 null 3.4209 and the sequence-blind pipeline 3.989 (S19). The system is doing real
   work; it is the gap between 3.21 and its own 2.31 A ceiling that the sprint is about.

There is NO published "ceiling for native-free peptide prediction at 9-16 residues" in the sense
the brief asks for -- no paper reports an upper bound, only method scores on small curated sets.
The nearest thing to a ceiling statement in the literature is McDonald 2023's negative result
(the best available confidence signal cannot rank within a peptide's own models), which is a
ceiling on SELECTION and agrees with the project's own measurement of the same thing.

## 5. Verdict table (topic 5)

| paper | information it contains that we do not | verdict |
|---|---|---|
| Maupetit et al. NAR 37:W498 (2009); Shen et al. JCTC 10:4745 (2014); Lamiable et al. NAR 44:W449 (2016); Rey et al. NAR 51:W432 (2023) -- PEP-FOLD 1-4 | the 27-state structural-alphabet profile as a RETRIEVAL KEY; sOPEP | KEPT only as the benchmark anchor (2.6 A selected at 9-25 aa). The retrieval key is REJECTED: closed here by S13's torsion-bin key (0.517 on FAIL18 vs a 0.562 majority baseline) and by the 22-key / 24-arm screen (no key beats BLOSUM on the chain; structural keys worsen pool-best and the ORACLE ceiling). sOPEP REJECTED (Legacy's functional class, measurably anti here) |
| Timmons & Hewage, Brief Bioinform 22:bbab308 (2021) -- APPTEST | a torsion head (closed here) and a fold-from-restraints operator (open, architectural) | KEPT as the 1.96 A anchor; REJECTED as a scorer |
| McDonald et al., Structure 31:111 (2023) | a controlled in-band ranking measurement at 10-40 aa | KEPT -- the field's only ceiling-like statement, and it is a ceiling on selection |
| Lindorff-Larsen et al., Science 334:517 (2011) | the equilibrium population as the selector | REJECTED on cost; KEPT as the statement of what the working selector is |
| Senior et al. Nature 577:706 (2020) -- AlphaFold1's distance potential | minimise the posterior directly rather than take a per-pair point estimate | KEPT as a route note (lanes M/X), not a method import |

## 6. What topic 5 tells the sprint

1. Do not quote 3.21 A against 1.96 A without the three caveats. The benchmarks differ in
   composition, the reporting differs (best-of-N vs single answer), and every published method
   selects within an ensemble its own energy generated.
2. Split the comparison. Generation: the project's 1.71-2.31 A ORACLE quantities are inside the
   published band, so generation at 9-16 aa is competitive. Selection: this is where the 0.9 to
   1.5 A sits, and the literature's own best method has no in-band skill here either.
3. The peptide literature has NO unexploited lead for this project. Its central architectural
   claim -- retrieve by predicted local conformation, not by sequence -- was tested here twice and
   failed both times, for a stated mechanism (section 3.1). I record this as a near-miss of my
   own: the S12 line that tempted me is a FAIL18 pool-best diagnostic, and the chain result that
   closes it is in a different file. Anyone reading S12's literature findings without the
   `docs/FINDINGS.md` screen beside it will re-propose this.
4. There is no published ceiling at this length. If the sprint produces one -- a measured bound
   on native-free selection for 9-16-mers with controls -- that is a contribution to the field
   and not merely to the project, because the field has not stated one.
