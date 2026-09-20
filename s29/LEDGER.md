# S29 LEDGER (append-only; `## S29-L<n> -- TITLE (date time, lane)`; number from the tail in one process)

## S29-L0 -- THE CHARTER, THE READING, THE POSITION, AND WAVE 1 (2026-09-19 23:35, coordinator)
Charter: `s29/BRIEF.md` (verbatim). Contract: `s29/S29_CONTRACT.md` (inherits S28's 12 rules and
its two addenda; adds rules 13 to 19: findings engaged per hypothesis, the information test for
every import, the nine questions and ten controls for quantum claims, 12-target probes before
126, the multiplicity count, mechanism beside outcome, the cost-RMSD meter). Running state:
`s29/STATE.md`. Governor v2.5 bands: suspend above 94% RAM, hard kill 95.5%, resume below 92%,
launch below 91% (the user's 94 to 95% instruction; `s26/governor.py` lines 57 to 65). Cap 6.

Reading done by the coordinator before this entry: `s27/REPORT_S28.md`, the S28 ledger tail,
`s27/RESUME_S28.md`, `s27/RETRACTIONS_S28.md`, lane A's damaged expectations,
`docs/STATE_BRIEF_2026-09-12.md` sections 4 and 5, `s26/REPORT.md` VII.2 (ESMFold infeasible on
this box, L13; the distogram already uses ESM-2 650M PCA-32, `core/data.py:886`), `s17`'s ESM
contact-head result (in-band content, length-gated, no shortlist value), `s23/LEDGER.md` L6 (the
scale is unreachable in principle), `s24/LEDGER.md` L2/L3 (score-selected pools have parallel
biases, cos 0.943; the unselected blind source is 31% independent and 0.76 A worse; mixtures sit
on a line). The "more information" routes the charter's second-order questions point at are
therefore already in the record: a larger language model (hardware), a structure-trained
predictor (leakage and hardware), the ESM contact head (measured), a learned residual generator
(S24), calibration of the posterior (S25 L2: worse under the readout), scale (S23 L6).

Position (`s29/STATE.md`): H1, the typicality axis (extrapolate the conditioned answer away
from the sequence-blind answer; the common-mode error as signal), with its strongest objection
already in the record (S24 L3's geometry says the sign is wrong) and a 20-minute ORACLE kill
(lane O rung 6). H0, the fallback that is probably the truth: the ceiling is the absence of
per-target information orthogonal to the typicality axis; S29 measures it decisively by
recognisable-operator class (lane O's ladder), searches for a source of such information (lanes
L, T, M), and builds the CVaR-VQE formulation in which the quantum stage consumes the pool
jointly so that the quantum stage is where any such information would act. H2, the divergent
direction (lane X): the pool is the wrong state space; a configuration-space CVaR-VQE with an
interacting Hamiltonian and a tail-ENSEMBLE readout, the one readout S13 to S21 never used.
Findings engaged: H1 attacks 3, 4, 5, 8, 11; H2 attacks 9, 10, 11 (aggregation) and accepts 8
as the risk; H0 accepts all eleven and measures 8 and 11.

Wave 1 (six lanes, spawned after this entry): L literature (permanent), M data-path map and
convenience choices and the harness audit, T theory (the contraction theorem, the locally
informative class, the spectral condition, CVaR over a structural observable, the reachable
set, the entropy term, what an objective must know), O the ORACLE ceiling ladder by operator
class and the typicality-axis probe, D adversary (permanent, rotating) and the cost-RMSD meter
and the suite, X the divergent state space. No build starts before O's ladder, T's sections 2
and 3, and M's map are posted. Budget: reading and theory 35%, probes 35%, full runs 30%.
Multiplicity: 0 endpoint comparisons at this entry.
Artefacts: this file; `s29/STATE.md`; `s29/briefs/S29{L,M,T,O,D,X}.md`.

## S29-L1 -- TOPIC 1, NATIVE-FREE STRUCTURE QUALITY: NO PUBLISHED QA METHOD HAS EVER BEEN TRAINED OR EVALUATED BELOW 40 TO 50 RESIDUES; THE FIELD'S FOUR SIGNAL CLASSES ARE EACH ABSENT, ALREADY OURS, LEAKAGE, OR MEASURED DEAD HERE; THE ONLY NATIVE-FREE SELECTORS SHOWN TO WORK AT 9 TO 25 aa RANK AN ENSEMBLE THEIR OWN ENERGY GENERATED (PEP-FOLD 2.6 A, APPTEST 1.96 A) OR ARE A CONVERGED FREE ENERGY (LINDORFF-LARSEN); AND AF2's pLDDT HAS NO WITHIN-TARGET SKILL ON 588 PEPTIDES OF 10 TO 40 aa, WHICH IS EXTERNAL CONFIRMATION OF CHARTER FINDING 8 (2026-09-19 23:54, L)

Question (brief topic 1, the coordinator's first need): is there ANY native-free quality signal
for a 9 to 16-mer that the field has shown to work, and on what data? Falsifier for the negative
half: a single published QA method trained or benchmarked at 9 to 16 residues, or any
general-purpose QA signal whose construction survives the length. No experiment; literature
only; 15 primary papers read for their construction, not their abstracts.

ANSWER. Two native-free selectors have peptide-length evidence, neither transfers.
(a) A coarse-grained force field ranking the ensemble IT generated: PEP-FOLD's lowest-sOPEP-
energy cluster centroid, 2.6 A mean cRMSD on 25 NMR peptides of 9 to 25 aa (Maupetit 2009);
APPTEST's XPLOR-NIH restraint energy, 1.96 A on 42 peptides of 9 to 25 aa (Timmons 2021).
Generation and selection share the energy, so the score is a convergence diagnostic over a
self-consistent set. Over a retrieval pool the same class of energy is measured WORSE than a
random subset here (AMBER +0.455, Legacy +0.330, 5/5 folds, S25 L16).
(b) The equilibrium population under converged all-atom MD (Lindorff-Larsen, Science 334:517,
2011; chignolin 10 aa, Trp-cage 20 aa), which is a FREE energy at 1e5 to 1e6 CPU-hours per
peptide; our budget is 6.43 core-equivalents and one AMBER process at a time.

THE NEGATIVE, WHICH IS THE POINT. Every general-purpose QA method excludes this length by
construction or by its own text: ProQ3 "targets that were shorter than 50 residues were filtered
out both from the CASP11 and CAMEO data sets"; VoroMQA's learning set is chains "longer than 99
residues" (12,825 entries); DeepAccNet trains on 50 to 300 residues; QMEANDisCo quotes its own
error bar as ~0.12 at 40 residues vs ~0.05 for larger; GraphQA and EnQA say nothing about short
chains. The four signal classes they reduce to (L_1 section 5): S1 packing/burial -- not a
VARIABLE at 9 to 16 aa (VoroMQA's own Table I separates good from bad by the solvent area share,
0.392 vs 0.470, which at peptide length is set by the length), and its CA forms ENV/HP/EXVOL are
coin tosses or anti on the built chain (S28-L48); S2 sequence-to-local-structure agreement -- IS
this project's axis (the distogram + DIS), anti at pref 0.071, and phi is sequence-blind at this
length (S13, 36.1 vs 36.4 deg); S3 consensus -- score_i = (1/N) sum_{j!=i} GDT-TS_ij
(DAVIS-EMAconsensus, verbatim from Kwon 2021) IS our CONS channel, and the family assumes
independent member errors while our pool is 68% common-mode by exact identity (S23 L9), with
CONS_POOL@chain the most anti-recognising of the 31 scorers (pref 0.056, pool-member contrast
-0.127, S28-L48); S4 single-point physics -- measured worse than random here.

THE EXTERNAL CONFIRMATION OF FINDING 8. McDonald et al., Structure 31:111-119 (2023), 588
peptides of 10 to 40 aa against NMR: "the lowest RMSD structures failed to correlate with lowest
pLDDT ranked structures"; full text, "there was no correlation between the first three ranks
assigned by AF2 and the structure that gave the lowest Ca RMSD" (preprint: 13% of lowest-RMSD
models were rank 0, ranks 1 to 4 took 26/18/24/19%). Taking rank-1 instead of best-of-5 costs
0.2 to 1.1 A, largest on the mixed-secondary-structure classes. That is the field's best
self-confidence signal, from the model that solved the protein case, with no within-target
ranking skill at peptide length -- the same failure this project measures, on a different
instrument, with vastly more information. Charter finding 8 is not a defect of the S27 library;
it is the field's position at this length.

ALSO RECORDED. (i) CASP14's assessors (Kwon 2021): single-model QA overtook consensus for the
first time (BAKER-experimental top-1 GDT-TS loss 8.4; BAKER-ROSETTASERVER top-1 lDDT loss 4.0;
DAVIS-EMAconsensus and GOAP "about average"), and "it is more difficult to train an EMA method
to estimate the superposition-dependent quantity, GDT-TS" -- our objective and our endpoint are
both superposition-dependent, which is the hard case by the field's own verdict. (ii) CASP15's
EMA winners are quasi-single-model (MULTICOM_qa, per-target corr 0.66, ranking loss 0.14): their
skill is bought with an INDEPENDENT reference ensemble, which this instrument does not have
(S24 L2/L3: the only second source is 31% angularly independent and 0.76 A worse). The value of
that family is the independence, not the scoring rule -- a statement about the generator, not
about QA. (iii) Training any of these leave-fold-out on our corpus is a dead end for three
measured reasons, not one: the features are already in `s27/cache`, a supervised combination of
all 31 reaches 0.960 held-out sign accuracy but the SAME rule on the RAND_SIGNED control reaches
0.952 with head-to-head 0.492 (it learned "is this a real trace", S28-L48), and in-band
discrimination here is signal-limited, not sample-limited (S13's flat learning curve).

THE ONE ACTIONABLE ITEM, and it is not an import. The only S4 sub-class with peptide-length
precedent is the free energy, and `docs/FINDINGS.md` section B records that S8 BUILT it --
`strain` = E(built) - E(freely relaxed), `E_free`, `F_qh` = <E> - kT S_qh with S_qh = 1/2 sum
log lambda of a 4 ps 300 K Langevin CA covariance, `F_boltz` = -kT log <exp(-E/kT)>, `width`,
`S_msf`, each through the same builder and relaxation as the native -- and it "completed only 1
of its 24 targets before the box filled up", committed and resumable as `python -m s8.relax
best`, with no claim made from one target. So the one class the literature says works at this
length is the one class this project started and never finished. Caveats stated with it: S_qh at
4 ps is noisy and downward-biased; the underlying channel (AMBER single point) is +0.455 A worse
than random, so the claim would be that the ENTROPY term rescues a channel measured anti; the
cost is one Langevin ensemble per candidate under the one-AMBER-process rule. Proposed as a
12-target pre-registered probe for lane D or O, NOT a build, and not by lane L (I run no
endpoint experiments).

VERDICT. Topic 1 yields no importable native-free scorer. KEPT (4, all as evidence or anchors):
Kwon 2021 (CASP14 EMA assessment); McDonald 2023 (the peptide pLDDT result); PEP-FOLD 2.6 A /
APPTEST 1.96 A as the published native-free anchors at 9 to 25 aa; Lindorff-Larsen 2011 (the
working selector is a free energy). REJECTED (5 families, with the reason each fails HERE):
packing/burial QA (the statistic is not a variable at this length); sequence-agreement QA (it is
our own axis, and it is anti); consensus QA (assumes independent member errors; 68% common-mode
by measurement); single-point physics QA (worse than random on this pool, with controls);
AF2-feature QA (leakage plus hardware, and its base signal fails at peptide length).
Multiplicity: 0 endpoint comparisons; no number in this entry is new measurement, every one is
cited to its sprint or its paper.
Artefacts: `s29/lit/L_1_native_free_qa.md` (15 papers with construction, assumption, check,
information test, verdict; sections 5 to 9 are the synthesis), `s29/lit/L_INDEX.md`.
