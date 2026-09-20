# L_1 -- NATIVE-FREE STRUCTURE QUALITY (charter finding 8; brief topic 1)

Lane L, 2026-09-20. The question the coordinator asked: is there ANY native-free quality signal
for a 9 to 16-mer that the field has shown to work, and on what data?

The instrument every assumption is checked against: 126 peptides of 9 to 16 residues; a
retrieval pool of 500 real CA windows per target with 68% common-mode error (S23 L9); a 17-bin
distogram from ESM-2 650M PCA-32 trained leave-fold-out on 6,800 windows; the built chain is the
ideal-geometry backbone rebuilt from phi/psi (`s12.instrument.project`); side chains only via
`sidechains.py` for AMBER; CPU only, 16 GB. The native-free library already measured on this
instrument: S27's channels (`s27/ham_lib.py`: RG_LAW, RG_UNIV, EXVOL, CAGEO, RAMA, CONTACT,
DISTPOT, ENV, HP, DSSPHB, ELEC, CONS, DMAP_CONS, TORS_CONS, POOLGO, SS_MATCH, CONTACT_LL,
DIS_MEAN, plus DIS, LEG and its 11 terms, AMB); S28-L48/L49: 31 scorers on the built chain, 20
prefer the projected production average to a 0.25 A ORACLE structure, 5 tie-dominated, 7 coin
tosses, the one two-clause pass (CONTACT@chain) sits at the max-over-31 null's mean.

The answer, stated first. No model-quality-assessment method in the literature has been trained
or benchmarked below about 40 to 50 residues; every one of them excludes or is silent on the
length regime of this instrument, and the signals they use fall into four classes, each of which
is either absent at 9 to 16 residues, already in the library in a stronger form, excluded by
the leakage rule, or measured dead here. At peptide length the only native-free signals ever
shown to select a structure are (i) a force field's lowest-energy conformer within an ensemble
the same predictor generated (PEP-FOLD's sOPEP, 2.6 A mean on 25 NMR peptides of 9 to 25 aa;
APPTEST's restraint energy, 1.96 A on 42 peptides of 9 to 25 aa) and (ii) the equilibrium
population under ms-scale MD (Lindorff-Larsen 2011), which is a free energy, not an energy.
AlphaFold2's pLDDT, the strongest self-consistency signal in the field, does NOT rank within its
own five models on 588 peptides of 10 to 40 aa (McDonald 2023). Details and verdicts below.

---

## 1. Single-model QA

### 1.1 VoroMQA -- Olechnovic K, Venclovas C. Proteins 85:1131-1145 (2017)

Construction (their eqs 1 to 10, read from the PDF). Atoms are typed (167 heavy-atom types),
contacts are Voronoi faces between atom cells, with solvent as an explicit partner type a0 and a
contact category set C (near/far by sequence separation x central/noncentral). The pseudo-energy
of a contact type is

    E(a_i, a_j, c_k) = log [ P_exp(a_i, a_j, c_k) / P_obs(a_i, a_j, c_k) ]              (1)

with P_obs the observed share of total contact AREA of that type in the learning set,
P_obs = S(a_i,a_j,c_k) / (S_int + S_sol), and P_exp the product of the marginals (a quasi-chemical
reference state: P_obs(a_i) P_obs(a_j) 2 P_obs(c_k) for i != j; P_obs(a_i) P_obs(c_0) for solvent).
An atom's normalised energy is the area-weighted mean over its contacts and its neighbours'
contacts, E_n(X) = sum E(type_x) area_x / sum area_x (eq 8); the atom score is
Q_a = 1/2 (1 + erf((E_n - mu_type)/(sigma_type sqrt 2))) (eq 9); the global score is the
BURIAL-WEIGHTED mean of atom scores, weight 1 for solvent-accessible atoms, 2 for their
neighbours, 3 for the next shell (eq 10).

Assumption. The discriminating statistic is packing: Table I of the paper gives the observed
probability of the SOLVENT category as 0.392 in the high-quality learning set against 0.470 in
CASP models, i.e. the good/bad separation is carried by how much area is buried. The learning
set is X-ray < 2.5 A, chains LONGER THAN 99 RESIDUES (12,825 entries). The paper states that
smaller structures receive lower global scores with greater variance (their Fig 3B), and gives
the usable rule as v < 0.3 bad, v > 0.4 good, the interval undecidable "using VoroMQA alone".

Check against this instrument. A 9 to 16-mer has no third burial shell and almost no second; the
burial weights collapse to 1, and the solvent share of a peptide's contact area is far above
0.47 whatever its fold. The statistic that separates good from bad in the paper's own table is
not a variable at this length: it is a constant set by the length. The project's CA-level
analogues (ENV = -log p(nbr bin | aa)/p(nbr bin), HP = -sum KD(aa) nbr, EXVOL) are the burial
statistic at the resolution this instrument has, and on the built chain they are coin tosses
(ENV@chain pref 0.492, HP@chain 0.456) or anti (EXVOL@chain 0.417 with 101 ties) (S28-L48).
VoroMQA also needs side chains; ours are rebuilt from the backbone by rotamer placement, so the
contact areas would be `sidechains.py`'s, not the candidate's.

Information test. Contact AREA with solvent as an explicit partner, instead of a distance
cut-off, and an atom-typed reference state. At peptide length the area of every contact is
dominated by the solvent face, and the atom typing rides on rebuilt side chains. Nothing
orthogonal to ENV/HP/EXVOL survives the length.

REJECTED: the signal it uses (buried area) does not exist as a variable at 9 to 16 residues; the
learning set starts at 100 residues.

### 1.2 ProQ2 / ProQ3 -- Uziela K, Wallner B, Elofsson A. Sci Rep 6:33509 (2016) (PMC5048106)

Construction. A per-residue regression to the S-score

    S_i = max(0, 1 - (d_i / d_0)^2),  d_0 = 3 A,                                         (2)

with the global score the mean of S_i over the chain. ProQ2 inputs: atom-atom and
residue-residue contact counts, surface area, agreement of the model's DSSP secondary structure
with PSIPRED, agreement of the model's RSA with predicted RSA, conservation from a PSI-BLAST
profile (three iterations against Uniref90), in windows of 5, 11 and 21 residues plus whole-model
values. ProQ3 replaces or adds Rosetta energies: full-atom (fa_atr, fa_rep, fa_sol, fa_elec,
hbond_sr_bb, hbond_lr_bb, hbond_bb_sc, hbond_sc, rama, omega, p_aa_pp, fa_dun, pro_close, ref)
and centroid (vdw, pair, env, cbeta, cenpack, rama; global rg, hs_pair, ss_pair, sheet, rsigma,
co). Training: CASP9, 117 targets, 30 models per target, average chain length 202.

Assumption. (a) A homology profile exists and predicted SS/RSA from it are informative; (b) the
Rosetta terms discriminate; (c) the regression transfers across targets. Test sets: "targets
that were shorter than 50 residues were filtered out both from the CASP11 and CAMEO data sets".

Check. (a) At this length a PSI-BLAST profile degenerates to the sequence; the project's
sequence channel (ESM-2 650M) is the stronger form of the same input, and S13 found phi
sequence-blind at peptide length (36.1 vs 36.4 deg). SS_MATCH (Chou-Fasman agreement) is the
library's form of the SS-agreement term: anti on the built chain (0.393). (b) The Rosetta
terms map onto LEG's 11 terms and AMBER, both measurably worse than a random subset (S25 L16).
(c) Chains under 50 residues are outside the method's own test domain.

Information test: none beyond the library; the profile-agreement idea is what the distogram is.
REJECTED (out of domain below 50 residues; every input class is already measured here).
### 1.3 ProQ4 -- Hurtado DM, Uziela K, Elofsson A. arXiv:1804.06281 (2018)

Construction (read from the PDF). Sequence branch: self-information I_i = -log(p_i / pbar_i)
and partial entropy S_i = -p_i log(p_i / pbar_i) per MSA column (their eq 1a/1b) plus one-hot
sequence, PRE-TRAINED on 5,687 PISCES chains to predict 3/6-state SS, RSA and sin/cos(phi,psi)
from sequence. Model branch: a COARSE description only: sin/cos(phi,psi), DSSP state (8 merged
to 6), RSA, DSSP backbone H-bond energies. Both branches enter a 1D fully-convolutional net;
fine-tuned on CASP9-10 models in a Siamese pair with a rank loss; target = per-residue lDDT,
global = mean. CASP11: R-local 0.77, R-per-model 0.56 (their Table I); per-target correlation
falls with mean model quality (their Fig 9, R 0.31).

Assumption. Quality is read as the agreement between what the sequence predicts locally
(SS, RSA, dihedrals) and what the model has locally. The paper's own reduced-input result says
this is where most of the signal is.

Check. This is exactly the axis the instrument already exhausts: the distogram is a
sequence-conditioned local-structure predictor and its agreement score (DIS) is the shipped
objective, which prefers the production average to the 0.25 A structure on 93% of targets
(S28-L48). A rank-loss learner over per-residue signed deviations is S13's set-transformer
(flat learning curve) and S28-L48's pairwise-logistic combination (anti-production, 0.492).

Information test: none. REJECTED (its signal is the project's own axis, already measured).

### 1.4 DeepAccNet -- Hiranuma N et al. Nat Commun 12:1340 (2021) (PMC7910447)

Construction. Predicts (i) per-residue Cb lDDT and (ii) the ESTOGRAM: for every residue pair,
a distribution of the SIGNED Cb-Cb distance error over 15 bins with edges -20,-15,-10,-4,-2,
-1,-0.5,0.5,1,2,4,10,15,20 A. Inputs: 3D convolutions over voxelised atoms in each residue's
local frame; 1D features (sequence, torsions, Rosetta intra-residue energies, SS); 2D features
(distances, orientations, Rosetta inter-residue energies); optional trRosetta distances and
ProtBert embeddings. Loss = estogram + 10 lDDT + 0.25 mask. Training: 7,307 X-ray chains,
50 to 300 residues, ~150 decoys each at GDT-TS 50 to 90 (RosettaCM, native perturbation,
trRosetta). Estograms are then used as pair restraints in Rosetta refinement (all-atom lDDT
+10% on 73 proteins). CASP14: the Baker single-model entries built on it were the best EMA
methods of the round (Kwon et al. 2021, section 3 below).

Assumption. A full-atom local environment (voxels, Rosetta energies) carries the error; decoys
are in the GDT-TS 50 to 90 band of a 50 to 300-residue chain; the paper notes the estogram
gets harder for larger proteins with more long-range interactions, and says nothing below 50.

Check. The voxel environment of a 9 to 16-mer is backbone plus rebuilt side chains and mostly
solvent; the Rosetta-energy inputs are the LEG/AMBER class that anti-ranks here. The estogram's
useful output, a native-free error DIRECTION per pair, was taken up in S16 (`s16/lit_FINDINGS.md`
B.1.2, E11) and every native-free error-direction arm chose "do nothing" (S16, closed).

Information test: the signed per-pair error distribution, which this project already tried in
its own basis. REJECTED (training domain 50 to 300 residues; the direction idea is closed by
S16 on this instrument).

### 1.5 QMEANDisCo -- Studer G et al. Bioinformatics 36:1765 (2020)

Construction. QMEAN single-model terms: torsion (phi/psi triplets), all-atom pair potential,
Cb pair potential, packing (atom count), agreement of DSSP SS with PSIPRED as a log-odds
S(d,p,c) = log[p(d,p,c) / (p(d) p(p,c))], agreement of accessibility with ACCpro. DisCo: from
HHblits homologues k, per pair (i,j) a Gaussian constraint g_ijk(d) = exp[-(d - mu_ijk)^2 / 2]
(the 1/2 is the published correction; width as printed) for template Ca-Ca distances under
15 A; homologues clustered, cluster weight w_c = exp[gamma SS_c], gamma = 70; s_ij(d) =
sum_c w_c h_ijc(d); DisCo_i = (1/n) sum_{j within 15 A} s_ij(d_ij). A feed-forward net weighs
DisCo against the QMEAN terms; target = per-residue lDDT; training CAMEO 9,500 models / 883
targets and CASP12 7,070 models / 70 targets. Per-residue AUC 0.87 (QMEAN) to 0.94 (with DisCo).
The paper states the expected error is ~0.12 for 40-residue models, converging to ~0.05 for
larger ones.

Assumption. The gain is the homologue restraint set; the single-model terms are the 2008 QMEAN
terms, which the paper needs DisCo to lift.

Check. Homologue distances are the leakage the fold clustering exists to exclude (identity-
clustered folds; `peptide_folds.json` pinned); at 9 to 16 residues HHblits hits are the fragment
bank the pool is already retrieved from. The single-model terms are RAMA, DISTPOT/CONTACT,
ENV, SS_MATCH here, all anti or coin tosses on the built chain (S28-L48). 40 residues is the
paper's own floor for a usable error bar.

Information test: template distances = the pool itself. REJECTED (leakage rule; below the
method's stated size floor; single-model half already measured).

### 1.6 GraphQA -- Baldassarre F et al. Bioinformatics 37:360 (2021)

Construction. A protein graph: nodes = residues with (one-hot aa, MSA self-information and
partial entropy as 23-dim vectors each, a 14-dim DSSP vector of dihedrals, RSA and SS type);
edges = pairs within d_max with features exp(-d^2/sigma) and a one-hot sequence separation in
classes {0,1,2,3,4,5:10,>10}. Message passing e'_ij = phi_e(e_ij, v_i, v_j, u); v'_i =
phi_v(rho(e'), v_i, u); targets local lDDT/CAD and global GDT_TS/TM/lDDT; loss = lambda_l l_l +
lambda_g l_g. Trained on 85k CASP9-12 decoys, ~270 targets; CASP13 R_target 0.779, local
Spearman 0.797, per-decoy Spearman 0.527.

Assumption and check. The inputs are the ProQ set in graph form; the ablation says DSSP matters
slightly more than MSA and that the raw-sequence variant is nearly as good (local RMSE 0.123 vs
0.121), i.e. the architecture is not where the signal is. Per-DECOY Spearman 0.527 is the
within-target number and it is the only one that matters here (S12: in-band is the only ranking
metric); it is obtained on CASP decoy sets spanning a wide quality range, not on an in-band
pool. No statement about short chains anywhere in the paper.

Information test: none (same inputs, better estimator). REJECTED.

### 1.7 EnQA -- Chen C et al. Bioinformatics 39:btad030 (2023)

Construction. SE(3)-equivariant GNN over the model's 3D graph, whose distinguishing input is
AlphaFold2's own output features (its predicted distogram/lDDT and structure) plus sequence and
structural features; target per-residue lDDT. Reported to beat AF2's self-reported pLDDT.

Assumption. An AF2 run on the target is available and its features are informative. That is a
structure-trained predictor on our targets, which is the project's leakage exclusion
(`docs/STATE_BRIEF` 5.7 / S26 VII.2: ESMFold is also infeasible on this box). Its own
advertised advantage rests on a quantity (pLDDT) that McDonald 2023 shows does not rank within
a peptide's own five models (section 4.1).

Information test: AF2's internal representation, which this instrument cannot have without
leakage. REJECTED (leakage + hardware; and the base signal fails at peptide length).

### 1.8 MD-based accuracy estimation -- Hu J, Xun S, Wu H, Wu Y, Jiang F, bioRxiv 439760 (2018)

Construction. Run short MD (< 1 ns) with the residue-specific force field RSFF2 and read
per-residue deviation/mobility as the accuracy estimate; 31 training and 24 test models;
claimed to reach single-model MAE state of the art, and to combine additively with a
knowledge-based score.

Assumption. Wrongly modelled regions are less stable under a good force field, i.e. the force
field ranks the native basin lower in free energy. On this instrument the force fields are
measured: AMBER as a ranker is +0.455 A WORSE than a random subset and Legacy +0.330 (S25 L16),
and AMBER-relaxation's apparent gain was against a worse projection (S16). The mechanism the
method needs (physics prefers the native) is the exact thing refuted here.
REJECTED (its assumption is falsified on this instrument, twice, with controls).

---

## 2. Consensus QA and its failure under correlated members

### 2.1 The naive consensus score (DAVIS-EMAconsensus; Pcons/ModFOLDclust lineage)

Construction, verbatim from the CASP14 assessment (Kwon S, Won J, Kryshtafovych A, Seok C,
Proteins 89:1940-1948, 2021): "DAVIS-EMAconsensus estimates model accuracy purely based on
consensus by scoring the ith model by an average GDT-TS to all other models in the pool as

    score_i = (1/N) sum_{j != i} (GDT-TS)_ij ."                                          (3)

Pcons is the same statistic with a different similarity; ModFOLDclust adds clustering.

Assumption. Errors are independent across models, so agreement is evidence of correctness. The
project has the exact-identity statement of how badly that fails here: the pool's error is 68%
common-mode (S23 L9, 50.7x the i.i.d. prediction), and "consensus is outlier avoidance, not a
nativeness signal" -- the native sits at the 82.8th percentile of the criterion the medoid
minimises (S12/S15). The literature states the same limit qualitatively: consensus "tends to
fail if the best models are far from the dominant structural cluster" (QMEANclust lineage,
recorded in `s15/LITERATURE.md` C16/N7), and Wang/Zhang's PWCom study (PLoS One 8:e74006, 2013)
gives CGDT = mean GDT to the rest of the set, notes it "works well when good models are among a
major cluster", and shows per-target cases where CGDT's Spearman is 0.39 while a single-model
potential reaches 0.89 (T0527) -- but with no quantitative threshold for when it flips.

Check. Eq (3) IS the project's CONS channel (mean CA-RMSD to the other members). On the built
chain CONS_POOL@chain has pref 0.056 and a pool-member contrast of -0.127: it is the single
most ANTI-recognising scorer in the 31 (S28-L48). Nine medoid/consensus arms closed 0% of the
in-pool gap (S12); the one positive in the record (score-filter + consensus medoid, -0.172 A,
S18) is outlier avoidance and is capped at the pool's mode.

Information test: none -- it is a within-pool statistic, and a bias shared by every member is
invisible to any within-pool statistic (the exact 68% identity). REJECTED, and this is the
family-level negative result the brief asks for: EVERY consensus QA method assumes decoy
errors are independent enough that agreement implies correctness; our pool violates the
assumption by construction, and the violation is measured, not argued.

### 2.2 Quasi-single-model QA (ModFOLD's -S variants, MULTICOM_qa, GATE)

Construction. Generate a reference ensemble for the target with an independent predictor, then
score the model by consensus AGAINST that ensemble rather than against its peers. In CASP15
this class (MULTICOM_qa, per-target correlation 0.66, ranking loss 0.14) topped the EMA table.

Assumption. The reference ensemble's errors are independent of the assessed model's. Here the
reference ensemble would have to come from the same retrieval pool and the same distogram, so
its errors are the same 68% common-mode bias: the contrast is zero by construction. S24 L2/L3
measured the nearest available approximation -- the unselected blind-library source is only 31%
angularly independent of the incumbent, and mixtures sit on a straight line (cos 0.943).

Information test: an independent predictor. That is the thing the instrument does not have and
cannot obtain native-free at this length (S26 VII.2 hardware; leakage for structure-trained
predictors). REJECTED as unavailable, and FLAGGED: if any lane finds a genuinely independent
generator, this is the QA family that would then apply -- its value is exactly the independence,
not the scoring rule.

---

## 3. What the CASP assessments actually say

Kwon S, Won J, Kryshtafovych A, Seok C. "Assessment of protein model structure accuracy
estimation in CASP14: old and new challenges." Proteins 89:1940-1948 (2021). Read from the
eScholarship manuscript.

- Classification rule: 70 methods, 46 single-model and 24 multi-model, separated by whether the
  score changes when other models are present (margin 0.02 on a 0 to 1 scale).
- The headline: "two single-model methods, BAKER-experimental and BAKER-ROSETTASERVER, performed
  better than the best multi-model method in this CASP, unlike CASP13 where the best multi-model
  method performed better than the best single-model method." Those are the DeepAccNet entries.
- Best top-1 losses: BAKER-experimental 8.4 GDT-TS units; BAKER-ROSETTASERVER 4.0 lDDT units.
  The reference methods DAVIS-EMAconsensus (eq 3) and GOAP (a statistical potential) were "about
  average".
- Absolute accuracy: best GDT-TS estimate was the NAIVE consensus (mean difference ~ 6.8); best
  lDDT estimate a single-model method (~6.7).
- "single-model methods estimate LDDT better than GDT-TS ... because it is more difficult to
  train an EMA method to estimate the superposition-dependent quantity, GDT-TS."
- "The difference in the top 1 GDT-TS loss between the top TS human method AlphaFold2 and the
  best EMA method is pronounced ... confirming that the current top tertiary prediction human
  group achieved results beyond possible from consensus."

Three things this instrument must take from it. (a) The field's best native-free quality signal
is a SUPERPOSITION-FREE local measure (lDDT); the superposition-dependent global measure, which
is what CA-RMSD is, is the one even the winners find hard -- our endpoint is the hard one.
(b) Every reported skill is measured on decoy sets that span a wide quality range; the project
already knows that garbage rejection reads as skill and does not transfer to in-band pools
(memory: decoy-bank-not-a-pool-proxy; in-band-is-the-only-ranking-metric). (c) CASP15's EMA
winners are quasi-single-model, i.e. they buy their skill with an independent predictor
(MULTICOM_qa: per-target correlation 0.66, ranking loss 0.14, first of 24 EMA predictors).

## 4. Peptide-specific: what works below 20 residues

### 4.1 AlphaFold2 on peptides -- McDonald EF, Jones T, Plate L, Meiler J, Gulsevin A, Structure 31:111-119 (2023); preprint bioRxiv 2022.02.17.480937

588 peptides, 10 to 40 aa, NMR references; classes: alpha-helical membrane / soluble, mixed
membrane / soluble, beta-hairpin, disulfide-rich. Published abstract, verbatim on the point
that matters: "AlphaFold2 showed several shortcomings in predicting Phi/Psi angles, disulfide
bond patterns, and the lowest RMSD structures failed to correlate with lowest pLDDT ranked
structures." The full text adds: "there was no correlation between the first three ranks
assigned by AF2 and the structure that gave the lowest Ca RMSD"; in the 155-peptide preprint
version only 13% of the lowest-RMSD structures were rank 0, with ranks 1 to 4 taking
26/18/24/19%. Best-of-5 CA-RMSD by class in the preprint: helical membrane 2.3 +/- 1.3, helical
soluble 4.5 +/- 2.8, mixed membrane 5.7 +/- 4.4, mixed soluble 4.4 +/- 2.8, beta-hairpin
2.9 +/- 1.4, disulfide-rich 2.2 +/- 1.2 A; the rank-1 model is 0.2 to 1.1 A WORSE than
best-of-5, the gap largest for the mixed classes. AF2 "performed better on short peptides than
long peptides".

Why this is the most important paper for finding 8. It is the field's largest controlled
measurement of a native-free ranker INSIDE a set of in-band peptide models, and the ranker is
the best self-confidence signal the field has, from the model that solved the protein case. It
has no within-target skill at 10 to 40 aa, and the price of its absence is 0.2 to 1.1 A -- the
same order as the gap this project is trying to close (3.21 vs the 2.31 A top-75 ORACLE
ceiling). This is independent, external confirmation that recognition, not generation, binds at
peptide length, and that it binds for a method with far more information than this instrument.
KEPT as evidence (not as a method): the external existence proof for S28's finding 8.

### 4.2 PEP-FOLD 1/2/3 -- Maupetit J, Derreumaux P, Tuffery P, NAR 37:W498 (2009); Shen Y et al. JCTC 10:4745 (2014); Lamiable A et al. NAR 44:W449 (2016)

Construction. Predict a profile over a 27-state structural alphabet of 4-residue fragments with
an SVM on a PSI-BLAST PSSM; assemble greedily under sOPEP (coarse-grained OPEP); cluster the
generated conformers (PEP-FOLD3: complete linkage on d = 1 - BCscore, cut at BCscore 0.8) and
rank clusters by sOPEP energy or Apollo.

Numbers. PEP-FOLD1 (2009): 25 NMR peptides, 9 to 25 aa; "averaged on 25 peptides and five runs,
the PEP-FOLD LEC reproduces the NMR structure at 2.6 A cRMSD" -- the LEC is the lowest-sOPEP-
energy cluster centroid, i.e. a DEPLOYABLE native-free selection. PEP-FOLD2 (2014): 56 peptides
of 25 to 52 aa, near-native for 95% of targets vs Rosetta 88%. PEP-FOLD3 (2016): near-native in
the top five scored models for 80% of targets. Neither the 2009 nor the 2016 paper reports a
correlation between sOPEP energy and RMSD, and the 2009 paper admits failures where the aqueous
experimental state is not the predicted one (magainin).

Assumption. sOPEP ranks within an ensemble that sOPEP ITSELF generated -- generation and
selection share the same energy, so the ensemble is already concentrated in sOPEP's low-energy
region and the selection is a local refinement of a self-consistent set. That is categorically
different from ranking a retrieval pool of real windows drawn from other proteins, where the
energy was never used to generate the members. S25 L16 is the matched measurement: as a ranker
over THIS pool a physics energy is worse than a random subset (AMBER +0.455, Legacy +0.330).

Information test. sOPEP contains a peptide-specific coarse-grained force field (backbone H-bond
cooperativity plus hydrophobic terms tuned on peptides) that this project does not have in
exactly that parameterisation -- but Legacy is the same functional class (11 terms including
hbond_local, hbond_longrange, coop_helix, coop_sheet) and is measured anti here.
REJECTED as a scorer for our pool (the ensemble-self-consistency assumption fails); KEPT as the
benchmark anchor for topic 5 (2.6 A native-free at 9 to 25 aa, 25 peptides).

### 4.3 APPTEST -- Timmons PB, Hewage CM, Brief Bioinform 22:bbab308 (2021)

Construction (recorded in `s12/lit_FINDINGS.md` Part 4): a 1D-CNN with gated residual blocks
predicts CA-CA/CB-CB distances AND phi/psi; folding by restrained simulated annealing in
XPLOR-NIH; final selection by lowest XPLOR-NIH energy / CYANA target function, deployable, not
oracle. Reported 1.96 A on 42 peptides of 9 to 25 aa (PEP-FOLD 2.05, PEPstrMOD 4.66).

Assumption for the selector. The ranked quantity is the RESTRAINT VIOLATION of a structure
against its own predicted distances -- again a self-consistency score over structures generated
to satisfy those same restraints. Our analogue is DIS (the L1 Bayes-risk score against the
distogram), which is the shipped objective and prefers the production average to a 0.25 A
structure on 93% of targets (S28-L48). The difference is not the score, it is that APPTEST's
candidates are generated BY the restraints (so the score is a convergence diagnostic) while ours
are retrieved independently (so the score is a discriminator, and it is a bad one).

Information test: a torsion head (phi/psi predicted directly) and a fold-from-restraints
operator with no retrieval bottleneck. The torsion half is closed here (S13: phi carries no
sequence signal at peptide length, 36.1 vs 36.4 deg blind; direct build emits 4.151 A). The
fold-from-restraints half is an open GENERATION route, not a QA route.
REJECTED as a QA signal; NOTED for lanes X/M as the one architecture at this length whose
terminal operator is not a pool average.

### 4.4 MD as the only true native-free selector at this length

Lindorff-Larsen K, Piana S, Dror RO, Shaw DE. Science 334:517-520 (2011). 12 fast-folding
proteins (chignolin 10 residues, Trp-cage 20), 100 us to 1 ms of all-atom MD with one force
field; the proteins "spontaneously and repeatedly fold to their experimentally determined native
structures". The selector is the equilibrium POPULATION -- a free energy that includes the
conformational entropy of the competing basins -- not a single-point potential energy.

Check. The project's record separates these cleanly: converged interaction-only AMBER flips the
sign on a decoy bank but does not transfer to real pools (memory:
amber-converged-interaction-only-ranks); single-point AMBER is +0.455 A worse than random here
(S25 L16). The quantity that works in Lindorff-Larsen is a Boltzmann average costing 10^5 to
10^6 CPU-hours per peptide; our budget is 6.43 core-equivalents and one AMBER process at a time.
REJECTED on cost by four to five orders of magnitude -- but it is the honest statement of what
the only demonstrated native-free peptide selector actually costs, and it says the missing
quantity is a FREE energy, which is a checkable structural statement about what any cheap
surrogate lacks.

---

## 5. The four signal classes, and why each is unavailable or dead at 9 to 16 residues

Every QA method above reduces to one or more of four signals. Stating them as classes is the
durable result, because it closes families rather than papers.

| class | what it measures | the method that uses it | status on this instrument |
|---|---|---|---|
| S1 packing / burial | how much area or how many neighbours are buried | VoroMQA, QMEAN packing, ProQ2 contacts, GOAP | NOT A VARIABLE at 9 to 16 residues (no second burial shell); its CA forms ENV/HP/EXVOL are coin tosses or anti (S28-L48) |
| S2 sequence-to-local-structure agreement | does the model's local conformation match what the sequence predicts | ProQ2/3/4, GraphQA, QMEAN SS/ACC agreement | IS the project's own axis: the distogram plus DIS. Anti on the built chain (DIS@chain pref 0.071), and phi is sequence-blind at this length (S13) |
| S3 consensus / agreement with other models | mean similarity to the rest of the pool | Pcons, DAVIS-EMAconsensus (eq 3), ModFOLDclust, QMEANclust | THE ASSUMPTION IS FALSIFIED HERE BY MEASUREMENT: 68% common-mode error (S23 L9); CONS_POOL@chain is the most anti-recognising of the 31 scorers (pref 0.056) |
| S4 physics energy or stability | force-field energy of the model, or its stability under MD | ProQ3's Rosetta terms, the MD-MAE method, sOPEP, XPLOR-NIH restraint energy | Measured worse than a random subset as a ranker on this pool (AMBER +0.455, Legacy +0.330, S25 L16); MD free energy works (Lindorff-Larsen) at 10^5 CPU-hours per peptide |

The only class with demonstrated peptide-length skill is S4, and only in its expensive
free-energy form, or in the degenerate form where the SAME energy generated the ensemble
(PEP-FOLD, APPTEST). There is a fifth class the field is now winning with -- an INDEPENDENT
predictor used as the reference ensemble (quasi-single-model QA, CASP15) -- and its whole value
is the independence, which this instrument does not have.

## 6. Could any of them be trained leave-fold-out on our corpus?

The brief asks this explicitly. Three obstacles, in order of severity.

1. The label. Every method above regresses to lDDT or the S-score of eq (2), computed against
   the native. Training leave-fold-out on our corpus is legitimate (the fold split is the same
   rule the distogram trains under), so the label is available for TRAINING windows. This is not
   the obstacle.
2. The features. All of S1 to S4 are already computed for every pool member in
   `s27/cache/<pdb>.npz`; a learner over them is exactly S28-L48's nested pairwise-logistic
   combination of 31 scorers. It reaches 0.960 held-out sign accuracy -- and the SAME rule on
   the RAND_SIGNED control reaches 0.952, head-to-head 0.492. It has learned "is this a real
   protein trace", not "is this one nearer the native". That is a measured, controlled negative
   for the whole supervised-QA-on-our-features programme.
3. The sample. S13's set-transformer over the full signed deviation map has a FLAT learning
   curve while a leaked label is loud at n=8 (memory: in-band-signal-limited-not-sample-limited).
   In-band discrimination here is signal-limited, not sample-limited, so more training data or a
   bigger QA architecture is the one thing that provably does not help.

Conclusion for the coordinator: importing any QA architecture from section 1 means training a
regressor on features the project already has, against a label the project can already compute,
in a regime where the last untested architecture class was measured flat. The expected value is
zero and the record says so three independent ways.

## 7. The answer to the coordinator's question

Is there ANY native-free quality signal for a 9 to 16-mer that the field has shown to work, and
on what data?

Yes, exactly two, and neither transfers:

(a) A coarse-grained force field ranking its OWN generated ensemble. PEP-FOLD's lowest-energy
    cluster centroid reaches 2.6 A mean cRMSD on 25 NMR peptides of 9 to 25 aa (2009), and
    APPTEST's XPLOR-NIH energy reaches 1.96 A on 42 peptides of 9 to 25 aa (2021). Data: small
    NMR peptide sets, 25 to 56 targets. It works because generation and selection share an
    energy; over a retrieval pool the same class of energy is worse than random here (S25 L16).

(b) The equilibrium population under converged all-atom MD (Lindorff-Larsen 2011, chignolin and
    Trp-cage at 10 and 20 residues). Data: 12 proteins, 100 us to 1 ms each. Cost is four to
    five orders of magnitude beyond this box.

And one strong negative from the largest peptide study in the field: AF2's pLDDT, the best
self-confidence signal available anywhere, has NO within-target ranking skill on 588 peptides of
10 to 40 aa (McDonald 2023); choosing its rank-1 model instead of its best costs 0.2 to 1.1 A.
No general-purpose QA method (VoroMQA, ProQ2/3/4, QMEAN(DisCo), GraphQA, DeepAccNet, EnQA) has
ever been trained or evaluated below 40 to 50 residues; ProQ3 filters out everything under 50
residues by name, VoroMQA's learning set starts at 100, DeepAccNet's at 50, QMEANDisCo quotes a
0.12 error bar at 40 residues as its small-model limit.

So: the project's finding 8 is not a local defect of its scorer library. It is the field's
position at this length, confirmed externally by the one method with far more information than
this instrument. A lane proposing to fix recognition by importing a QA method is proposing
something no published method supports at 9 to 16 residues.

## 8. What this leaves open (for the coordinator, not a recommendation)

Three things in the literature are NOT closed by the above, and each is a different lane:

1. Independence, not scoring. CASP15's EMA winners are quasi-single-model: their skill is bought
   with a reference ensemble whose errors are independent of the assessed model's. The project's
   measured obstacle is that its only second source is 31% angularly independent and 0.76 A
   worse (S24 L2/L3). The literature says the VALUE of an independent generator is large; the
   project says the one it has is too correlated. That is a statement about the generator, not
   about QA.
2. lDDT, not RMSD. The CASP assessors state plainly that superposition-free local measures are
   the learnable ones and GDT-TS (hence CA-RMSD) is the hard one. This project's objective and
   endpoint are both superposition-dependent. Whether a superposition-free local objective
   (a per-pair or per-window lDDT-like quantity) behaves differently in the cost-RMSD meter is
   an unasked question here, and it is cheap for lane D to ask.
3. The free-energy gap, and it is NOT unexplored -- it is UNDER-POWERED and declared as such.
   The only working peptide selector in the literature is a free energy (4.4), and every physics
   channel the project has SHIPPED is a single point. But `docs/FINDINGS.md` section B records
   that S8 built exactly the right stage -- `strain` = E(built) - E(freely relaxed); `E_free`
   (energy at the candidate's own relaxed minimum); `F_qh` = <E> - kT S_qh with
   S_qh = 1/2 sum log lambda of the superposed CA covariance of a 4 ps 300 K Langevin ensemble;
   `F_boltz` = -kT log <exp(-E/kT)>; plus `width` and `S_msf` -- computed through the same
   builder and relaxation as the native, seeded from the restrained-relaxed candidate so it
   samples that candidate's own basin. It "completed only 1 of its 24 targets before the box
   filled up", it is committed and resumable (`python -m s8.relax best`), and the record
   explicitly says no free-energy claim is made from one target. The one recorded note is
   directional: on 1A13 the TOTAL energy puts the native at the 74th percentile while the
   INTERACTION-ONLY energy puts it at the 28th.
   So the single S4 sub-class with peptide-length precedent in the literature is the one class
   this project started, never finished, and never closed. Caveats before anyone gets excited:
   the S_qh estimate is noisy and downward-biased at 4 ps; the ranker underneath is AMBER, which
   is +0.455 A worse than random as a single point (S25 L16), so the claim would have to be
   that the ENTROPY term rescues a channel measured anti, which is a strong claim; and the cost
   is one Langevin ensemble per candidate under the one-AMBER-process rule. It belongs to lane D
   (cost meter) or O (ORACLE ladder) as a 12-target probe, not to a build.

## 9. Verdict table (topic 1)

| # | paper | signal class | information it contains that the project does not | verdict |
|---|---|---|---|---|
| 1.1 | VoroMQA, Olechnovic & Venclovas, Proteins 85:1131 (2017) | S1 | contact AREA with explicit solvent partner; atom-typed quasi-chemical reference | REJECTED -- burial is not a variable at 9-16 aa; learning set >99 residues |
| 1.2 | ProQ2/ProQ3, Uziela et al., Sci Rep 6:33509 (2016) | S1+S2+S4 | none (PSSM profile is weaker than ESM-2; Rosetta terms = LEG/AMB class) | REJECTED -- targets <50 residues filtered out by the authors |
| 1.3 | ProQ4, Hurtado et al., arXiv:1804.06281 (2018) | S2 | a Siamese RANK loss on lDDT over a coarse model description | REJECTED -- its axis is our DIS axis; rank-learners over our features are measured (S13, S28-L48) |
| 1.4 | DeepAccNet, Hiranuma et al., Nat Commun 12:1340 (2021) | S1+S2+S4 | the estogram: a signed per-pair error DISTRIBUTION, usable as restraints | REJECTED -- trained 50-300 residues; the error-direction idea is closed here by S16 |
| 1.5 | QMEANDisCo, Studer et al., Bioinformatics 36:1765 (2020) | S1+S2+S3 | homologue-derived per-pair distance constraints with identity weighting | REJECTED -- homologue distances are the leakage rule and, at this length, the pool itself; authors' own floor ~40 residues |
| 1.6 | GraphQA, Baldassarre et al., Bioinformatics 37:360 (2021) | S1+S2 | none (same inputs, graph estimator; ablation says architecture is not the signal) | REJECTED |
| 1.7 | EnQA, Chen et al., Bioinformatics 39:btad030 (2023) | S2 (AF2 features) | AF2's internal representation | REJECTED -- leakage + hardware; and its base signal (pLDDT) fails at peptide length (4.1) |
| 1.8 | MD-MAE, Hu et al., bioRxiv 439760 (2018) | S4 | short-MD stability as an absolute per-residue accuracy | REJECTED -- assumes physics prefers the native, refuted here with controls (S25 L16, S16) |
| 2.1 | DAVIS-EMAconsensus / Pcons / ModFOLDclust; PWCom, Wang & Zhang, PLoS One 8:e74006 (2013) | S3 | none -- it is our CONS channel, eq (3) | REJECTED -- family-level: consensus QA assumes independent member errors; our pool is 68% common-mode by measurement |
| 2.2 | quasi-single-model QA (ModFOLD-S, MULTICOM_qa, GATE) | S3 against an independent ensemble | an INDEPENDENT predictor as the reference | REJECTED as unavailable; FLAGGED -- the value is the independence, not the scoring rule |
| 3 | Kwon et al., CASP14 EMA assessment, Proteins 89:1940 (2021) | -- | the field's own verdict: superposition-free local measures are learnable, GDT-TS/RMSD is not | KEPT as evidence (method-free) |
| 4.1 | McDonald et al., Structure 31:111 (2023) | -- | a controlled in-band ranking measurement at 10-40 aa with the best confidence signal in the field | KEPT as evidence -- external confirmation of finding 8 |
| 4.2 | PEP-FOLD 1/2/3 | S4 (self-consistent) | sOPEP, a peptide-tuned coarse-grained force field | REJECTED as a pool scorer (self-consistency assumption); KEPT as the 2.6 A native-free anchor at 9-25 aa |
| 4.3 | APPTEST, Timmons & Hewage, Brief Bioinform 22:bbab308 (2021) | S4 (self-consistent) | a torsion head + a fold-from-restraints terminal operator | REJECTED as QA; NOTED as a generation architecture (1.96 A at 9-25 aa) |
| 4.4 | Lindorff-Larsen et al., Science 334:517 (2011) | S4 (free energy) | the equilibrium POPULATION, i.e. an entropy term | REJECTED on cost (10^5-10^6 CPU-h/peptide); KEPT as the statement of what the working selector actually is |

KEPT (4, all as evidence or anchors, none as an importable operator): Kwon 2021; McDonald 2023;
PEP-FOLD's 2.6 A / APPTEST's 1.96 A anchors; Lindorff-Larsen 2011's free-energy statement.
REJECTED families (5): packing/burial QA; sequence-agreement QA; consensus QA; single-point
physics QA; AF2-feature QA.

Nothing in topic 1 is importable as a new native-free scorer for this instrument. The one
actionable item is section 8.3 (the unfinished S8 free-energy stage), and it is a probe for
another lane, not a literature import.
