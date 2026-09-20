# S29 LANE L -- LITERATURE INDEX (one line per paper: verdict, and the information it adds)

Format: topic file :: citation :: signal class :: what it contains that the project does not ::
KEPT / REJECTED. Signal classes are defined in `L_1_native_free_qa.md` section 5
(S1 packing/burial, S2 sequence-to-local-structure agreement, S3 consensus, S4 physics/stability).
The instrument every check is made against: 126 peptides of 9-16 residues, a 500-member
retrieval pool with 68% common-mode error, an ESM-2 650M PCA-32 17-bin distogram trained
leave-fold-out, built chain as the reporting basis, CPU only, 16 GB.

## Topic 1 -- native-free structure quality (`L_1_native_free_qa.md`, ledger S29-L1)

| paper | class | information it adds | verdict |
|---|---|---|---|
| Olechnovic K, Venclovas C. VoroMQA. Proteins 85:1131-1145 (2017) | S1 | Voronoi contact AREA with solvent as an explicit partner type; atom-typed quasi-chemical reference | REJECTED -- burial is not a variable at 9-16 aa; learning set is chains >99 residues |
| Uziela K, Wallner B, Elofsson A. ProQ3. Sci Rep 6:33509 (2016) | S1+S2+S4 | none (PSSM < ESM-2; Rosetta terms = the LEG/AMB class) | REJECTED -- the authors filter out every target under 50 residues |
| Hurtado DM, Uziela K, Elofsson A. ProQ4. arXiv:1804.06281 (2018) | S2 | a Siamese RANK loss on per-residue lDDT over a coarse model description | REJECTED -- same axis as DIS; rank-learners over our features are measured flat (S13, S28-L48) |
| Hiranuma N et al. DeepAccNet. Nat Commun 12:1340 (2021) | S1+S2+S4 | the estogram: a signed per-pair error distribution usable as a restraint | REJECTED -- trained on 50-300 residues; error-direction closed here by S16 |
| Studer G et al. QMEANDisCo. Bioinformatics 36:1765 (2020) | S1+S2+S3 | homologue-derived per-pair distance constraints, identity-weighted (gamma=70) | REJECTED -- leakage rule; at this length the homologues ARE the pool; authors' floor ~40 residues |
| Baldassarre F et al. GraphQA. Bioinformatics 37:360 (2021) | S1+S2 | none (graph estimator over the ProQ feature set; ablation says architecture is not the signal) | REJECTED |
| Chen C et al. EnQA. Bioinformatics 39:btad030 (2023) | S2 | AlphaFold2's internal features as QA input | REJECTED -- leakage + hardware; base signal fails at peptide length |
| Hu J et al. MD-based MAE. bioRxiv 439760 (2018) | S4 | sub-ns MD stability as an absolute per-residue accuracy | REJECTED -- assumes physics prefers the native; refuted here with controls (S25 L16) |
| Kwon S, Won J, Kryshtafovych A, Seok C. CASP14 EMA assessment. Proteins 89:1940-1948 (2021) | -- | the field's own verdict: single-model QA overtook consensus in CASP14; lDDT is learnable, GDT-TS (hence RMSD) is not | KEPT as evidence |
| Wang Q, Zhang Y et al. PWCom. PLoS One 8:e74006 (2013); Pcons/ModFOLDclust/DAVIS-EMAconsensus lineage | S3 | none -- score_i = mean similarity to the rest of the pool IS our CONS channel | REJECTED (family-level: consensus assumes independent member errors; our pool is 68% common-mode) |
| quasi-single-model QA (ModFOLD-S, MULTICOM_qa, GATE; CASP15) | S3 | an INDEPENDENT predictor as the reference ensemble | REJECTED as unavailable; FLAGGED -- the value is the independence, not the scoring rule |
| McDonald EF et al. Benchmarking AlphaFold2 on peptide structure prediction. Structure 31:111-119 (2023) | -- | a controlled in-band ranking measurement on 588 peptides of 10-40 aa: pLDDT has NO within-target skill; rank-1 costs 0.2-1.1 A vs best-of-5 | KEPT as evidence -- external confirmation of charter finding 8 |
| Maupetit J, Derreumaux P, Tuffery P. PEP-FOLD. NAR 37:W498 (2009); Shen Y et al. JCTC 10:4745 (2014); Lamiable A et al. NAR 44:W449 (2016) | S4 (self-consistent) | sOPEP, a peptide-tuned coarse-grained force field; 27-state structural-alphabet profile retrieval | REJECTED as a pool scorer (it ranks its OWN ensemble); KEPT as the 2.6 A native-free anchor at 9-25 aa |
| Timmons PB, Hewage CM. APPTEST. Brief Bioinform 22:bbab308 (2021) | S4 (self-consistent) | a torsion head and a fold-from-restraints terminal operator (no retrieval bottleneck) | REJECTED as QA; NOTED as a generation architecture (1.96 A at 9-25 aa) |
| Lindorff-Larsen K, Piana S, Dror RO, Shaw DE. Science 334:517-520 (2011) | S4 (free energy) | the equilibrium population: an entropy term, not a single-point energy | REJECTED on cost (10^5-10^6 CPU-h per peptide); KEPT as the statement of what the only working peptide selector is |

Topic 1 running count: 15 entries, 4 KEPT (all as evidence or anchors; none importable as an
operator), 11 REJECTED across 5 families.

## Topic 2 -- correlated error in ensembles (`L_2_correlated_error.md`, ledger S29-L8)

| paper | what it needs that we lack | verdict |
|---|---|---|
| Krogh A, Vedelsby J (1995), ambiguity decomposition, as Brown, Wyatt & Tino, JMLR 6:1621-1650 (2005) eq (10) | nothing -- it is an identity, and it IS S23 L9's decomposition (160.36 + 63.82, exact to 2.7e-14) | KEPT as framing: finding 11 is a law, not a defect |
| Ueda N, Nakano R (1996), bias-variance-covariance, same paper eq (9): bias^2 + (1/M)var + (1-1/M)covar | nothing | KEPT: the (1-1/M) coefficient caps aggregation; yields the M-table and the derived <= ~0.008 A bound on "more members" |
| Liu Y, Yao X (1999) negative correlation learning; bounds lambda_upper = M/(M-1), gamma_upper = M^2/(2(M-1)^2) (Brown et al. eq 39) | members that are ESTIMATORS BEING TRAINED | REJECTED -- ours are retrieved real windows; noted for any future parameterised generator (lane X) |
| Abe T, Buchanan EK, Pleiss G, Cunningham J. Pathologies of Predictive Diversity in Deep Ensembles. arXiv:2302.00704 | nothing; it is a negative result (R_ens = R_avg - Jensen gap) | KEPT as confirmation -- diversity interventions harm good ensembles; matches "diversity-maximising selection: dead" |
| Glynn PW, Szechtman R. Some New Perspectives on the Method of Control Variates. MCQMC 2000, Springer (2002) | a control variate whose mean E Y is KNOWN | REJECTED -- the known mean is the native; the formal reason within-pool statistics cannot see the common mode |
| Peherstorfer B, Willcox K, Gunzburger M. SIAM Review 60:550-591 (2018), eqs (3.12), (3.16) | m0 samples of the HIGH-FIDELITY quantity as an unbiased anchor | REJECTED -- it reduces variance around an anchor; it never de-biases a biased model |
| Richardson-style extrapolation / two-source contrast (classical; H1's family) | a KNOWN bias ratio along a SHARED direction | REJECTED as stated -- S24 L2/L3 gives cos 0.647 with the second source 0.76 A worse; the non-parallel 31% is amplified |
| Boosting / gradient residual fitting | a residual learnable from inference features, with decorrelated errors | REJECTED -- closed by S24; `error-coherence-decides-correctors` gives the violated condition (+0.31 A coherent vs -0.14 A i.i.d. at equal accuracy) |
| Jumper J et al. AlphaFold2. Nature 596:583-589 (2021) -- recycling | training the predictor with its own output in the loop | REJECTED (training-time, leakage, hardware); KEPT as evidence: pLDDT fits lDDT at r = 0.76 globally yet has no in-band skill at peptide length |

Topic 2 running count: 9 entries, 4 KEPT (framing/confirmation), 5 REJECTED families.
Derived in the note (arithmetic on `s23/results/errdecomp.json`, not a new measurement):
an infinite pool of the same kind returns 3.040 A against the shipped 3.0483 point cloud, so
"average more members" is worth <= ~0.008 A.

## Topic 3 -- decision theory of structure point estimates (`L_3_decision_theory.md`, ledger S29-L12)

| paper | what it adds | verdict |
|---|---|---|
| Classical Bayes estimators (L2 -> posterior mean, L1 -> median, 0-1 -> mode), Berger ch. 4; with the Jensen contraction \|\|E X - E Y\|\| <= E\|\|X - Y\|\| | the one-line proof that a coordinate average must shrink every interatomic distance | KEPT as framing: the measured 25.8% contraction is Jensen, not a distogram bias, and cannot be reweighted away |
| Blau Y, Michaeli T. The Perception-Distortion Tradeoff. CVPR 2018 (arXiv:1711.06077), Thm 1 and Thm 3 | for ANY distortion measure, the distortion-optimal estimator's output distribution must diverge from the real one, most steeply at low distortion | KEPT -- the framing theorem for charter finding 8: realism-type scorers must disprefer the RMSD-optimal answer |
| the same paper's converse (perfect realism by drawing unrelated real signals) | realism alone buys no accuracy | KEPT -- the theoretical form of the pool-member control (S28-L36/L37) |
| Guo C, Pleiss G, Sun Y, Weinberger KQ. On Calibration of Modern Neural Networks. ICML 2017 (arXiv:1706.04599) | temperature scaling is argmax-invariant (and median-invariant for symmetric posteriors) | REJECTED as a lever; KEPT as the mechanism behind S25 L2's "calibrating makes RMSD worse" |
| Gneiting T, Raftery AE. JASA 102(477):359-378 (2007) | properness (the distogram's cross-entropy IS the strictly proper log score); CRPS as the distance-sensitive diagnostic | REJECTED as a training change; KEPT as the correct posterior-comparison metric, in place of MAE |
| Senior AW et al. Nature 577:706 (2020); Jumper J et al. Nature 596:583 (2021) | the field's answer to an over-confident per-pair posterior was to stop taking a per-pair point estimate | KEPT as a route note supporting lane X; not importable |

Topic 3 running count: 6 entries, 5 KEPT (framing/diagnostic/route), 2 levers REJECTED.
Answer to the brief's specific question ("is the Bayes estimator of an over-confident posterior
contracted?"): NO, not from over-confidence per se -- a symmetric width error moves neither the
mean nor the median. The contraction here is averaging (Jensen) plus shrinkage toward the prior's
centre (the typicality axis).

## Topic 4 -- quantum (`L_4_quantum.md`, ledger S29-L13)

| paper | what it adds | verdict |
|---|---|---|
| Barkoutsos PKl, Nannicini G, Robert A, Tavernelli I, Woerner S. Quantum 4:256 (2020), arXiv:1907.04769, eqs (11) and (12) | CVaR is DEFINED on sorted samples and is scoped to "classical optimization problems, which yield diagonal Hamiltonians" | KEPT -- the definitional source of the project's set-equality theorem; cite it rather than presenting S25's theorem as a discovery |
| the same paper: the overlap result and Prop 5.1 | CVaR's global-minimiser set is {theta : overlap >= alpha}, large and flat; estimator SE grows as 1/alpha; CVaR and expectation landscapes do not map to each other | KEPT and FLAGGED -- a candidate mechanism for "the objective is optimised but the structure does not move", and a new falsifier for our own quantum arms |
| Cerezo M, Sone A, Volkoff T, Cincio L, Coles PJ. Nat Commun 12:1791 (2021), arXiv:2001.00550 | locality of the OBSERVABLE, not depth alone, governs gradient decay (assumes local 2-designs) | NOTED, not actionable -- our ansatz is nowhere near a 2-design; contract rule 9 |
| Larocca M et al. Quantum 6:824 (2022), arXiv:2105.14377 | the dynamical Lie algebra as the diagnostic of reachability and gradient scaling | KEPT and RECOMMENDED -- the project's DLA is listed as never measured; classical linear algebra, minutes, answers charter question 5 |
| Ragone M et al. Nat Commun 15 (2024), arXiv:2309.09342 | an exact loss-variance expression through the DLA; variance falls as dim(g) grows | KEPT with a sourcing caveat (abstract read; equation not quoted verbatim) |
| Cerezo M et al. Nat Commun 16:7907 (2025), arXiv:2312.09121 | trainable, barren-plateau-free models tend to admit classical simulation of the loss | KEPT -- the framing constraint: "it trains" is now evidence FOR simulability; the charter's ten controls are this paper's programme |
| Amin MH, Andriyash E, Rolfe J, Kulchytskyy B, Melko R. Phys Rev X 8:021050 (2018) | a thermal state of a non-commuting H whose measurement distribution is not a classical Boltzmann distribution | KEPT -- the cleanest instance of (C1) non-commuting AND (C2) not-an-eigenvector together |
| QAOA's mixer (structural reading) | an off-diagonal term is the only thing that moves amplitude BETWEEN candidates | KEPT -- explains S28's rank-one failure as a degenerate mixer, not as "off-diagonals do not help" |
| Grimsley HR, Economou SE, Barnes E, Mayhall NJ. ADAPT-VQE, Nat Commun 10:3007 (2019), arXiv:1812.11173 | grow the ansatz one operator at a time by the commutator gradient over an operator pool | KEPT conditionally -- the principled way to build a circuit for a structured non-diagonal H; not itself a source of non-classicality |

Topic 4 running count: 9 entries, 7 KEPT, 1 NOTED, no RMSD lever. The condition for
non-classicality: (C1) the Hamiltonian's terms do not commute AND (C2) the prepared object is not
an eigenvector (Gibbs/thermal, or a sampling distribution, or a free-energy objective). The
project has never satisfied both.

## Topic 5 -- peptide prediction at 9-16 residues (`L_5_peptide_ceiling.md`, ledger S29-L14)

| paper | what it adds | verdict |
|---|---|---|
| Maupetit J, Derreumaux P, Tuffery P. NAR 37:W498 (2009) -- PEP-FOLD1 | 2.6 A cRMSD on 25 NMR peptides of 9-25 aa for the lowest-sOPEP-energy cluster centroid (a genuine native-free single answer) | KEPT as the published anchor |
| Timmons PB, Hewage CM. Brief Bioinform 22:bbab308 (2021) -- APPTEST | 1.96 A on 42 peptides of 9-25 aa, selected by XPLOR-NIH energy; plus a fold-from-restraints terminal operator | KEPT as the anchor and as an architectural note (lanes M/X); REJECTED as a scorer |
| Lamiable A et al. NAR 44:W449 (2016); Shen Y et al. JCTC 10:4745 (2014); Rey J et al. NAR 51:W432 (2023) -- PEP-FOLD2/3/4 | best-of-5 hit rates at 25-52 aa; PEP-FOLD4 ties trRosetta and AF2 on 17 peptides | KEPT as context; the retrieval key REJECTED (closed here, see below) |
| McDonald EF et al. Structure 31:111 (2023) | the field's only ceiling-like statement, and it is a ceiling on SELECTION | KEPT |
| Lindorff-Larsen K et al. Science 334:517 (2011) | the equilibrium population as the working selector at 10-20 aa | REJECTED on cost; KEPT as the statement of what the working selector is |
| Senior AW et al. Nature 577:706 (2020) -- AlphaFold1 | minimise a potential fitted to the distance posterior instead of taking a per-pair point estimate | KEPT as a route note |

Topic 5 running count: 6 entries. NO published ceiling exists at this length. The comparison
splits: our GENERATION (ORACLE 1.71-2.31 A) is inside the published band; our SELECTION is the
0.9-1.5 A gap. **Closed, and nearly re-imported by this lane**: PEP-FOLD's central architectural
claim (retrieve by predicted local conformation) failed here twice -- S13's torsion-bin key
(0.517 on FAIL18 vs a 0.562 majority baseline) and the 22-key / 24-arm screen (no key beats
BLOSUM on the chain; structural keys worsen pool-best 1.711 -> 2.161 and the ORACLE ceiling
1.994 -> 2.542 by concentrating the pool on the predicted type).

## Topic 5 -- peptide prediction at 9-16 residues (`L_5_peptide_ceiling.md`)
(pending)
