# Sprint 15 — LITERATURE, COMPARABILITY MATRIX, AND NOVELTY BOUNDARY

Agent: `lit`. Date: 2026-09-05. Extends `s14/lit_FINDINGS.md` (which located and read
arXiv:2609.02113 / QTF in full). **This document does not repeat that audit; it extends it.**
Where a fact was established in Sprint 14 from a primary source, it is cited as `[S14]` and not
re-fetched. New primary fetches this sprint are marked `[P15]`.

Method: WebSearch / WebFetch only. No repository compute. No benchmark artefact touched.

**Tiering** per BRIEF §5: **LITERATURE-SUPPORTED** (primary source fetched and read this sprint
or in S14), **LITERATURE-REPORTED (UNVERIFIED)** (search snippet / abstract / secondary only),
**HYPOTHESIS** (my inference, not stated by the source), **REFUTED**.

**Fabrication discipline:** every citation below carries a verification mark. Nothing is
asserted as a number unless I read it in a fetched source or it is carried forward from S14's
own verification list. Where I could not confirm, the row says UNVERIFIED and says what
specifically is missing.

---

## STATUS: COMPLETE. Sections 1-6 written.

Sections completed so far are marked DONE. This file is written incrementally so it survives
interruption.

---

## 1. HEADLINE FINDINGS OF THIS SPRINT (read this first)

Four results change the Sprint 14 novelty picture. Two of them take claims we held.

### 1.1 CLAIM-TAKING: the certified optimum being worse than random IS PUBLISHED — on a lattice contact potential

> Roget M, Damour C, Cadet F, Wang J. **"Assessing Cost Hamiltonian Reliability in Quantum
> Protein Structure Prediction."** arXiv:2606.21241v1 [quant-ph], 19 June 2026.
> Affiliations: University Paris City & University of Reunion; Centre for Quantum Information,
> Simulation and Algorithms, University of Western Australia; ENERGYLab, University of Reunion;
> PEACCEL. Code: `github.com/mroget/Paper-Code--Assessing-Cost-Hamiltonian-Reliability-in-Quantum-Protein-Structure-Prediction`

**[P15] — I fetched the PDF and read the complete body, Methods, Results and Conclusion.**

This paper does, on a tetrahedral lattice with a contact potential, the thing Sprint 14
believed nobody else could do: **it enumerates the entire feasible space and reports the RMSD
of the cost function's certified global minimum, and finds it worse than a random draw.**

Verbatim from their Conclusion (LITERATURE-SUPPORTED):

> "On average, the error of the structure with minimal cost is larger than the error of a
> structure taken at random. Thus, for short peptides, the minimum-cost conformation has, on
> average, a larger RMSD than a randomly selected feasible conformation. This result does not
> indicate a failure of quantum optimization algorithms per se. Rather, it shows that, for the
> cost Hamiltonian studied here, the ideal solution of the encoded optimization problem may be
> structurally inaccurate. **The bottleneck is therefore not necessarily the quantum optimizer,
> but the objective function encoded into the Hamiltonian.**"

And from Results §4.1 (LITERATURE-SUPPORTED):

> "The minimum RMSD quantifies the degradation in accuracy due to using a tetrahedral lattice.
> This degradation is not negligible at around 1.5 A. However, we see that the RMSD of the
> solution with minimum cost is much higher at around 5 A... Even worse, the average RMSD
> across all solutions is, on average across all peptides, better than the optimal solution
> according to the cost Hamiltonian."

**Their protocol, exactly (all LITERATURE-SUPPORTED):**

* Dataset: **12,446 peptides (4,866 unique sequences), length 6-600 aa, from the PDB.** Exact
  sequence duplicates removed; one representative per sequence chosen deterministically as the
  first `(pdbid, chain)` in lexicographic order. Entries with missing residues or missing CA
  were discarded. They repeated the analysis with duplicates and saw the same qualitative trends.
* Representation: **self-avoiding walks on a tetrahedral lattice**, arc length 3.8 A, `n-1`
  turns from `{1,2,3,4}`. CA only; all-atom reconstruction via `pulchra` is downstream and
  not part of the analysis.
* Cost function: shell-resolved contact energy,
  `p(c,s) = sum_{i<j} sum_{m=1..l} 1[|i-j| > k_min and (m-1)d_max < ||c_i - c_j||_2 <= m d_max] * lambda_m * eps_{s_i,s_j}`
  with `d_max = 7.8`, `k_min = 5`, `lambda = (1)` for the first-order study. `eps` = Miyazawa-Jernigan.
* **Enumeration: complete, for every peptide of length <= 15.** For length > 15 they sample
  10,000 feasible solutions with the Pivot algorithm; the 95% CI amplitude on Spearman rho at
  10,000 samples is **approximately 0.05** (their Appendix A).
* Metric: CA-RMSD after alignment, against the experimental structure. Per peptide they report
  average RMSD over the space, minimum RMSD over the space, RMSD of the minimum-cost solution,
  and Spearman rho(cost, RMSD).
* Correlation study restricted to length >= 13, "because the solution spaces of these peptides
  contain too many solutions with the same cost" below that — **i.e. they hit the tie problem
  our own BRIEF flags (tie-averaged argmin) and solved it by exclusion rather than by
  tie-averaging.**

**Their length result, which is new information for us and which we should adopt:**

> "One can see that the correlation is negative for small peptides, null for small proteins of
> length 50 and becomes positive only for larger proteins. This is not surprising since the MJ
> energy coefficients used in the cost function were calculated using large protein (of length
> higher than 50)."

**Our instrument (9-16mers) sits exactly in the regime where they measure rho(cost, RMSD) as
NEGATIVE.** LITERATURE-SUPPORTED. This is the first external, quantitative statement that the
energy-ranking failure is a *short-peptide* pathology that *attenuates with chain length*.
It reframes our negative results: they may be a property of the length regime rather than of
energy functions in general.

**Sensitivity analysis they ran, which matters for our Legacy-vs-AMBER work:** replacing MJ
contact energies with an HP-type model gave "comparable" correlations and the same qualitative
trend with interaction shells. Their reading: "once the problem is projected onto a coarse CA
lattice representation, the structure of the encoded objective -- including contact
definitions, interaction shells, and sequence-separation constraints -- may dominate over the
fine details of the pairwise coefficient matrix." **This is the external analogue of our
"Legacy is 98.8% one term / a clash gate" finding.**

**What this takes from us, and what it does not.** See §5. In one line: *"the certified global
optimum of the encoded objective is worse than random"* is **TAKEN for a coarse-grained
contact potential on a lattice** and is **NOT taken for an all-atom force field (AMBER ff14SB)
or for an off-lattice torsion representation**, which is what we measure. Our claim must be
re-scoped, not withdrawn — and we should now **cite Roget et al. as corroboration** rather than
present the phenomenon as unobserved.

### 1.2 CORROBORATION: decoy-based scoring-function evaluation is known to be artefact-ridden

> Handl J, Knowles J, Lovell SC. **"Artefacts and biases affecting the evaluation of scoring
> functions on decoy sets for protein structure prediction."** *Bioinformatics* 25(10):1271-1279
> (2009). doi:10.1093/bioinformatics/btp150

**[P15] — abstract/full-text page fetched; the quoted sentences were returned verbatim.**

* Verbatim: **"For 139 out of 149 of the decoy sets considered, such trivial discrimination is
  indeed possible, as the native corresponds to an extremum under at least one Amber99 term."**
  (~93% of decoy sets.)
* They separate three bias sources verbatim: **"(i) a lack of independence between sampling
  points; (ii) the use of the knowledge of the native during decoy generation; and (iii) the
  use of scoring functions during decoy generation."**
* Verbatim: **"evaluation based on the rank/z-score of the native is a weak test of scoring
  function performance"**, and correlation analysis can be "significantly affected by biases
  intrinsic to particular decoy sets."

**This is a direct external confirmation of our memory entry "Decoy-bank ranking skill does not
transfer to real pools"** (garbage rejection reads as skill), published seventeen years ago.
LITERATURE-SUPPORTED. We should cite Handl et al. 2009 whenever we defend the decision to
confirm every ranker on real pools rather than decoy banks. It also means our decoy-bank
finding is **not novel** and should be presented as a replication, not a discovery.

### 1.3 CORROBORATION AT THE STATE OF THE ART: AlphaFold2 cannot rank its own five peptide models

> Gulsevin A, Meiler J. **"Benchmarking Peptide Structure Prediction with AlphaFold2."**
> bioRxiv 2022.02.17.480937 (preprint v1). Published as **"Benchmarking AlphaFold2 on peptide
> structure prediction", *Structure* 31(1) (Jan 2023), PMID 36525975.**

**[P15] — bioRxiv v1 full text fetched.**

From the **preprint v1** (LITERATURE-SUPPORTED): 155 peptides, 16-60 aa, PDB entries with
solid-state or solution NMR structures. Metric: CA-RMSD after CA alignment, against the
lowest-energy reference NMR model. Their accuracy bands: <3 A "good", 3-6 A "moderate",
>6 A "poor". Mean +- sd CA-RMSD by class:

| class | mean RMSD (A) | without outliers |
|---|---|---|
| helical membrane-associated | 2.3 +- 1.3 | 2.2 +- 1.0 |
| helical soluble | 4.5 +- 2.8 | 2.4 +- 0.9 |
| mixed SS, membrane | 5.7 +- 4.4 | 3.5 +- 0.8 |
| mixed SS, soluble | 4.4 +- 2.8 | 2.9 +- 1.1 |
| beta-hairpin | 2.9 +- 1.4 | 2.7 +- 1.2 |
| disulfide-rich | 2.2 +- 1.2 | 2.1 +- 1.0 |

**The selection result, verbatim and load-bearing:**

> "Only 13% of the lowest-RMSD structures came from structures that had a rank of 0 (best
> rank). For the remaining predictions, 26%, 18%, 24%, and 19% came from ranks 1, 2, 3, and 4
> respectively."

and their conclusion:

> "our results suggest that the pLDDT metric... is not a meaningful metric to classify peptide
> conformations generated by AF2."

**Read that against a uniform null.** Five models; uniform selection would put 20% of the
lowest-RMSD structures at each rank. They measure 13 / 26 / 18 / 24 / 19. **AlphaFold2's own
confidence ranking of its own five peptide models is at or slightly below chance.** This is the
state of the art in structure prediction failing exactly the way our objectives fail, on
peptides, published in *Structure*. LITERATURE-SUPPORTED (the percentages are theirs; the
comparison to a 20% uniform null is my arithmetic, and trivial).

**Version discrepancy, flagged, not resolved.** A search snippet for the published *Structure*
version describes **588 peptides of 10-40 aa** and "length normalized Ca RMSD", where the
preprint I fetched describes 155 peptides of 16-60 aa and plain CA-RMSD. I fetched the
preprint; ScienceDirect returned HTTP 403 for the journal version. **The 588 / 10-40 /
length-normalised figures are LITERATURE-REPORTED (UNVERIFIED)** and must be re-checked before
citing. The 155 / 16-60 figures and both verbatim quotations are LITERATURE-SUPPORTED.

### 1.4 SCOPE-NARROWING: the closest predictive comparator to our instrument is 4.89 A

> Zhang Y, Yang Y, Cheng F, Lu C-C, Saeidi N, Volchenboum SL, Zhao J, Chen S, Jiang W, Guan Q.
> **"A Hybrid Quantum-AI Framework for Protein Structure Prediction on NISQ Devices."**
> arXiv:2510.06413v1, 7 October 2025.

**[P15] — HTML full text fetched.**

* Dataset: **75 protein fragments, 10-14 residues, from PDBbind.** 5 VQE candidates each =
  **375 conformations.**
* Quantum part: VQE on IBM 127-qubit hardware; Hamiltonian `H_geom + H_steric + H_int + H_misc`
  with **Miyazawa-Jernigan** pair energies; CA-only backbone candidates.
* Classical part: **NetSurfP-3.0** supplies SS3/SS8 probabilities, predicted phi/psi and
  relative solvent accessibility as statistical potentials. (**Correction to the S14 ledger,
  row 13**, which wrote "NSP3": the paper's object is NetSurfP-3.0. Same model, better name.)
* Selection: `E_fuse(c) = alpha*Ehat_q(c) + beta*Dhat_ss(c) + gamma*Dhat_angle(c)`, candidates
  sorted ascending, **"The top-ranked structure c(1) is selected as the final prediction."**
  **This is a genuinely predictive, native-free selection rule.** LITERATURE-SUPPORTED.
* Result: **mean RMSD 4.89 A, median 4.70 A, sd 1.10 A** across the 75 fragments.
* **No CVaR.** Not mentioned anywhere. (Confirms the S14 correction to the Sprint 13 ledger.)
* RMSD definition: "typically computed over backbone or Ca atoms" — **atom set and chain range
  are not pinned down.** UNVERIFIED.

**Why this is the closest comparator we have, and why it is still not a fair one.**
10-14 residues overlaps our 9-16mers; the selection is native-free; the metric is CA-ish. But
their pool is **5 candidates** where ours is 75-500; their targets are PDBbind binding-site
fragments, a biased structural class; their RMSD atom set and chain range are undefined; and
there is no fold-disjointness or identity-clustering statement anywhere. **Our 3.204 A
(`synthesis_fit`, 126 targets, full-chain Kabsch CA, no oracle) is numerically better than
their 4.89 A, but the comparison is indicative only and must carry all four caveats.** See §4.

---

## 2. NEGATIVE EVIDENCE — who else has seen what Sprint 14 saw

The BRIEF asked for this to be as strong as the positive literature. It is.

### 2.1 "The objective optimum is not the structural optimum" — well attested, across four settings

| # | source | setting | the negative, in their words or numbers | ver |
|---|---|---|---|---|
| N1 | **Roget, Damour, Cadet, Wang. arXiv:2606.21241 (Jun 2026)** | tetrahedral lattice, MJ contact energy, 12,446 peptides 6-600 aa, **complete enumeration at length <= 15** | min-cost RMSD ~5 A vs lattice-best ~1.5 A; **"the minimum-cost conformation has, on average, a larger RMSD than a randomly selected feasible conformation"**; rho(cost,RMSD) **negative for small peptides**, ~0 at length 50, positive only for larger proteins | **[P15]** |
| N2 | **Cumbo et al. arXiv:2609.02113 (QTF, Sep 2026)** | off-lattice all-atom torsions, 7.3 M structures, 3 energy backends | **"OpenMM shows the opposite tendency: total OpenMM potential energy is negatively correlated with RMSD in both proteins"**; native-like hit rate falls from 99-100% in the pool to **0.3-6% after energy selection**; "model selection remains the main bottleneck"; **no numeric rho printed anywhere** | **[S14]** |
| N3 | **Gulsevin & Meiler, bioRxiv 2022 / *Structure* 2023** | AlphaFold2, peptides | only **13%** of lowest-RMSD structures were AF2's own rank-0 model, against a **20% uniform null**; "pLDDT... is not a meaningful metric to classify peptide conformations generated by AF2" | **[P15]** |
| N4 | **Handl, Knowles, Lovell. *Bioinformatics* 25:1271 (2009)** | 149 decoy sets, Amber99 | **"For 139 out of 149 of the decoy sets considered, such trivial discrimination is indeed possible, as the native corresponds to an extremum under at least one Amber99 term"**; "evaluation based on the rank/z-score of the native is a weak test of scoring function performance" | **[P15]** |
| N5 | **Maffucci & Contini, *JCTC* 12(2):714-727 (2016)** (see §2.4) | REMD, **six AMBER ff (ff96, ff99SB/ildn/ildn-phi series — NOT ff14SB) x three implicit solvent models**, 8 peptides | verbatim: **"a combination of the force field and implicit solvation models able to accurately predict the native structure of all the considered peptides was not identified"** | **[P15]** abstract |
| N6 | coarse-grained energy-landscape literature | CG protein models | "non-native structures are observed with energies similar to those of the native state"; CG models struggle particularly "for small peptides" | **[S]** |
| N7 | CASP model-quality-assessment literature (QMEANclust lineage) | consensus/clustering QA | consensus methods "tend to fail if the best models are far from the dominant structural cluster" | **[S]** |

**N7 is the external statement of our memory entry "Consensus is outlier avoidance, not a
nativeness signal — the mechanism is capped at the pool mode."** LITERATURE-REPORTED
(UNVERIFIED — snippet only). Worth a primary fetch before we cite it in a paper.

**N1 is the most important row in this document.** See §5 for what it takes.

### 2.2 "Optimisation loses to random sampling" — one strong precedent, and it is on our problem

> Boulebnane S, Lucas X, Meyder A, Adaszewski S, Montanaro A. **"Peptide conformational
> sampling using the quantum approximate optimization algorithm."** *npj Quantum Information*
> **9**, 70 (2023). arXiv:2204.01821.

**[P15] — arXiv abstract page fetched, abstract quoted verbatim below. The full-text
quantification of "a small overhead" was NOT obtained. That number is UNVERIFIED.**

> "...based on numerical simulations on 20 qubits, we find less promising results: deep quantum
> circuits are required to achieve accurate results, and **the performance of QAOA can be
> matched by random sampling up to a small overhead.** Overall, these results cast serious
> doubt on the ability of QAOA to address the protein folding problem in the near term, even in
> an extremely simplified setting."

**This is the closest published relative of Sprint 14 Result 1** ("running VQE is worse than
not running it, against best-of-the-same-number-of-shots from the untrained circuit"). The four
differences all matter:

* Their control is **uniform random sampling of conformations**. Ours is **best-of-N shots from
  the SAME untrained parameterised circuit** — a strictly tighter control, because it holds the
  ansatz induced distribution fixed and isolates the *optimisation step* rather than the
  *circuit*. **I found no published instance of the untrained-circuit best-of-N control.** §5.
* Their algorithm is QAOA; ours is CVaR-VQE over a torsion library.
* Their result is *parity* ("matched... up to a small overhead"). Ours is **a loss**: 0/12,
  +0.65 to +1.32 A, including from ORACLE warm starts. A loss is a stronger statement than a tie.
* Their setting is an "extremely simplified" lattice model; ours is off-lattice with a real
  force field.

**Do not overclaim.** A hostile referee will say "Boulebnane already showed QAOA is matched by
random sampling." The correct answer is that matched-by-random is not the same measurement as
loses-to-its-own-untrained-circuit, and only the second identifies *concentration onto a bad
objective tail* as the mechanism. **The mechanism is the surviving claim; the parity
observation is not.**

### 2.3 Oracle-selected reporting — how common? ANSWER: endemic, and almost nobody names it

The BRIEF asked specifically for this. **LITERATURE-SUPPORTED as a pattern** (the table below is
built only from sources I or S14 fetched); **the meta-claim that it is unnamed is HYPOTHESIS**,
supported by a targeted and unsuccessful search for a paper that names, quantifies or
standardises the practice.

| paper | oracle-selected headline | predictive number | is the gap named? |
|---|---|---|---|
| **QTF, arXiv:2609.02113** | **0.623 A** chignolin, argmin **by RMSD to native** over 2.68 M snapshots; 2.501 A Trp-cage | **2.90 A / 5.39 A** median final model | Both printed and correctly labelled in Table 3, but **the abstract leads with the oracle number**, and the gap is never discussed as a methodological problem. [S14] |
| **Kannan et al. arXiv:2510.15316** | **"The minimum RMSD structures predicted using quantum computers consistently lie in the range 1.22 A to 3.11 A"** — the operative word is *minimum* | none reported | **No.** The only accuracy statement in the abstract is a minimum over predicted structures. **[P15]** abstract |
| **Roget et al. arXiv:2606.21241** | reports min-RMSD over the enumerated space, but **explicitly as a representational ceiling**: "the solution closest from ground truth that can be expressed on the tetrahedral lattice... measures the degradation of the prediction due to using a tetrahedral lattice" | RMSD of the min-cost solution, labelled **"the best result one can hope to obtain by running the quantum optimization"** | **YES. This is the one paper that gets it right**, and it separates representational ceiling from selection outcome by construction. **[P15]** |
| **Zhang et al. arXiv:2510.06413** | none — selection is by fused energy | **4.89 A mean** | N/A, no oracle used. **[P15]** |
| **Gulsevin & Meiler** | the lowest-RMSD-of-5 analysis | rank-0 (pLDDT-selected) performance | **YES, implicitly** — the 13%-at-rank-0 result *is* the gap, and they draw the right conclusion from it. **[P15]** |

**Three conclusions for our positioning.**

1. **Kannan et al. (arXiv:2510.15316) is the number our "< 2.0 A" goal will be read against, and
   its 1.22 A is a minimum over predicted structures on a lattice, not a predictive mean.**
   Anyone placing it beside our full-chain, no-oracle, 126-target mean is comparing
   incomparable quantities. The matrix in §4 says so explicitly, per the BRIEF.
2. **Roget et al. independently arrived at the discipline we use** — separate "best achievable
   in this representation" from "what the objective returns" — three months before this sprint.
   **Cite them for the protocol; do not claim it.**
3. The gap **is not named as a phenomenon in any source I could find.** A short methodological
   contribution that (i) names oracle-selected reporting, (ii) tabulates the gap across the QPSP
   literature, and (iii) proposes `pool_best / selected / returned` as a mandatory reporting
   triple, appears **genuinely open** — and it is nearly free for us, because
   `s12/instrument.py` already emits exactly that triple
   (`pool_best 1.711 / top75_best 2.306 / synthesis_fit 3.204`). HYPOTHESIS, well supported.

### 2.4 Force-field failure modes on peptides — ~~RESOLVED AT ABSTRACT LEVEL, AND IT PARTLY LETS US OFF~~ **THE DOWNGRADE IS REVERSED (Sprint 16)**

> ## ⚠ SPRINT 16 CORRECTION — 2026-09-06, RETRACT workstream. THE ff14SB THREAT IS REINSTATED.
>
> **The abstract of the *PCCP* 2018 study flagged below as "UNVERIFIED" and "the highest-value
> unfetched source in this document" has now been obtained verbatim**, via the Europe PMC REST API
> (`europepmc.org/webservices/rest/search?query=DOI:"10.1039/c7cp08010g"&resultType=core`) after RSC
> returned 403 and PubMed a cookie wall. Tier: **[A+] — the publisher-deposited abstract record, not
> a search snippet and not a summarising model's rendering.** Full text still not obtained
> (subscription; `isOpenAccess: N`, `inEPMC: N`), so Methods-level detail below is inference and is
> labelled as such.
>
> > Shao Q, Zhu W. **"Assessing AMBER force fields for protein folding in an implicit solvent."**
> > *Phys. Chem. Chem. Phys.* **20**(10):7206–7216 (2018). doi:10.1039/c7cp08010g. PMID 29480910.
> > 39 citations.
>
> **Verbatim, from the abstract record:** *"we performed enhanced sampling MD simulations to assess
> the ability of six AMBER force fields (FF99SBildn, FF99SBnmr, FF12SB, FF14ipq, **FF14SB**, and
> **FF14SBonlysc**) as coupled with a recently improved pair-wise **GB-Neck2** model in modeling the
> folding of two helical and two β-sheet peptides. Whilst most of the tested force fields can yield
> roughly similar features for equilibrium conformational ensembles and detailed folding free-energy
> profiles for short α-helical **TC10b** in an implicit solvent, the measured counterparts are
> significantly discrepant in the cases of larger or β-structured peptides (**HP35, 1E0Q, and GTT**).
> Additionally, the calculated folding/unfolding thermodynamic quantities can only partially match
> the experimental data. Although a combination of the force fields and GB-Neck2 implicit model able
> to describe all aspects of the folding transitions towards the native structures of all the
> considered peptides was not identified, we found that **FF14SBonlysc coupled with the GB-Neck2
> model seems to be a reasonably balanced combination** to predict peptide folding preferences."*
>
> **This is our exact force-field / implicit-solvent pair.** The condition on which §2.4 rested —
> *"this paper does NOT test ff14SB"* — is true of Maffucci & Contini and **false of Shao & Zhu**.
>
> ### What was actually tested, and the length regime (systems pinned against the PDB)
>
> | system | residues | class | verdict in the abstract |
> |---|---|---|---|
> | **TC10b** (Trp-cage variant, cf. 2JOF) | **20** | α-helical | force fields agree; *this one works* |
> | **HP35** (villin headpiece, cf. 1YRF) | **35** | α-helical (3-helix) | "significantly discrepant" |
> | **1E0Q** (N-terminal 17-mer of ubiquitin) | **17** | β-hairpin | "significantly discrepant" |
> | **GTT** (FiP35 variant, Pin1 WW domain) | ~35 | β-sheet | "significantly discrepant" |
>
> Residue counts for 1E0Q, 2JOF and 1YRF taken from the RCSB entry API
> (`deposited_polymer_monomer_count`); the TC10b↔2JOF and GTT↔FiP35 identifications are **[I]**.
>
> ### THE RE-SCOPING, AND WHETHER IT IS DEFENSIBLE — the hostile reading first
>
> LIT proposes re-scoping on the grounds that Shao & Zhu measure *folding thermodynamics from
> enhanced-sampling MD*, not *conformer ranking of a fixed pool*. **That defence is real but weaker
> than it looks, and two of its natural supports fail outright.**
>
> **What fails.**
> * **"They only tested proteins longer than ours" — FALSE.** `1E0Q` is a **17-residue** β-hairpin,
>   one residue above our 9–16 band, and it is one of the systems on which the force fields are
>   reported discrepant. The length defence is not available for the β-sheet failure.
> * **"Free energy is a different, weaker observable than potential-energy ranking" — BACKWARDS.**
>   Their observable is *stronger*: a folding free-energy profile is Boltzmann-weighted and includes
>   the entropy our single-point and minimised-pool rankings ignore. If ff14SB/GB-Neck2 does not put
>   the native basin at the free-energy minimum, a potential-energy ranking of a pool drawn from the
>   same landscape has no reason to do better. **Their result is upstream of ours, and it predicts
>   ours.**
>
> **What holds, and is worth stating.**
> * Their conclusion is **weaker than "ff14SB fails"**: no *single* combination described *all
>   aspects* of *all four* peptides, the α-helical 20-mer is reported as working, and
>   **FF14SBonlysc/GB-Neck2 is endorsed as "a reasonably balanced combination"**. A referee quoting
>   this paper against us is quoting a mixed result, not a clean refutation. *(LIT's paraphrase —
>   "finds no combination describes folding to the native" — over-reads it, and is corrected here.)*
> * They report **no conformer ranking**, no enumerated space, no certified global optimum, no
>   per-target CI, and no filtering/ordering decomposition. Nothing in the paper says where a native
>   sits in the *rank order* of a fixed candidate pool.
> * They test **4 peptides**; we test 126 with pinned folds, and 19 with a certified optimum over
>   1.28e7 exactly labelled structures.
>
> ### THE CORRECTED POSITION
>
> **Cite Shao & Zhu as establishing the expectation; claim only the certification.** The honest form
> is: *"ff14SB with GB-Neck2 is already documented not to describe folding to the native for β-sheet
> peptides down to 17 residues (Shao & Zhu 2018). We do not claim to discover that AMBER mis-ranks
> peptide natives — we quantify it in a different observable (rank position within a fixed pool),
> over an enumerated space with a certified global optimum, with paired fold-aware CIs on 126
> targets, and we extend it to polarizable physics (AMOEBA) which does not fix it."*
>
> **This project's own record is consistent with Shao & Zhu on every point, which is the strongest
> reason to treat the re-scoping as honest rather than self-serving:**
> * native anti-ranking is **distributed across terms**, not a single-term defect — i.e. exactly the
>   "no combination fixes all aspects" shape;
> * **AMOEBA**, a polarizable force field, nearly flips trpzip but fails chignolin on vdW — the
>   ceiling extends past fixed-charge parameterisation, so it is not an ff14SB-specific bug;
> * **converged interaction-only AMBER flips the sign on the decoy bank but does not transfer to
>   real pools** — garbage rejection, not nativeness;
> * Sprint 15 PHYS: the **accuracy** claim (−0.022 Å) is **absent on the frame-reproducible 69% of
>   the instrument** (−0.009 [−0.024, +0.006], 38W/46L, ×3 frames), while the **validity** claim is
>   CONFIRMED and larger than recorded (Ramachandran 0.466 → **0.874**, 116W/2L; clashes 1.397 →
>   **0.000**, 63W/0L).
>
> **Net effect on the novelty ledger.** "AMBER anti-ranks peptide natives" moves from *finding* to
> *confirmation of documented expected behaviour*. What survives is (i) the **observable** — rank
> position in a fixed pool against a certified optimum, which Shao & Zhu do not measure; (ii) the
> **statistics**; (iii) the **AMOEBA extension**; and (iv) the **positive** result, which was never
> threatened by this literature at all: AMBER's demonstrated role here is **stereochemical repair**,
> and that half is large, replicated and unclaimed by anyone.
>
> **STILL NOT OBTAINED:** the *PCCP* full text. It remains the highest-value fetch for Methods-level
> questions (which RMSD/Q measure, which sampling protocol, whether any *ranking* statistic appears).
> Its **abstract-level** conclusion is now settled and no longer needs fetching to be quoted.

**CORRECTION, same session.** An earlier draft of this section said this source was unfetched
and that it covered ff14SB. **Both were wrong. Correction retained per BRIEF §1.5.**

> Maffucci I, Contini A. **"An Updated Test of AMBER Force Fields and Implicit Solvent Models
> in Predicting the Secondary Structure of Helical, beta-Hairpin, and Intrinsically Disordered
> Peptides."** *J. Chem. Theory Comput.* **12**(2):714-727 (2016).
> doi:10.1021/acs.jctc.5b01211. PMID 26784558.

**[P15] — abstract obtained and quoted verbatim. Full text NOT obtained (PubMed cookie wall,
ACS and RSC both HTTP 403).**

Verbatim abstract: *"Replica exchange molecular dynamics simulations were performed to test the
ability of six AMBER force fields and three implicit solvent models of predicting the native
conformation of two helical peptides, three beta-hairpins, and three intrinsically disordered
peptides."*

Verbatim conclusion: *"**Although a combination of the force field and implicit solvation models
able to accurately predict the native structure of all the considered peptides was not
identified**, the GB-Neck2 model seems to well compensate for some of the conformational biases
showed by ff96 and ff99SB/ildn/ildn-phi. The force fields of the ff99SB series coupled with
GB-Neck2 reasonably discriminated helices from disordered peptides, while a good prediction of
beta-hairpin conformations was only achieved by performing two independent simulations: one
with the ff96/GB-Neck2 combination and the other with GB-Neck2 coupled with any of the
ff99SB/ildn/ildn-phi force fields."*

**THE IMPORTANT SCOPE POINT: this paper does NOT test ff14SB.** The named force fields are
**ff96 and the ff99SB / ff99SB-ildn / ff99SB-ildn-phi series**, with GB-Neck2 among the solvent
models. So it does **not** directly document our exact `core/amber.py` pair (ff14SB + GBn2).

**Consequences, both directions:**

* ~~**It does NOT close our claim 3 the way I feared in an earlier draft.** The documented
  failure is for ff96/ff99SB-series peptide *secondary-structure* prediction, not for ff14SB
  *conformer ranking*. Our AMBER anti-ranking result is therefore **not** simply a restatement
  of a known ff14SB defect. **The threat to claim 3 from this source is downgraded.**~~
  **⚠ THE DOWNGRADE IS REVERSED, Sprint 16, 2026-09-06 (RETRACT).** The statement remains true
  *of Maffucci & Contini*, but the sibling source below (Shao & Zhu, *PCCP* 2018) **does** test
  ff14SB and FF14SBonlysc with GB-Neck2 — our exact pair — on peptides including a **17-residue**
  β-hairpin, and reports the folding thermodynamics "significantly discrepant". See the correction
  box at the head of this section. The threat is **reinstated and re-scoped, not withdrawn**.
* **It still constrains the framing.** "AMBER-family force fields with GB implicit solvent do
  not reliably place the native lowest for peptides" is established for adjacent members of the
  family, and a referee will say so. **Cite Maffucci & Contini in our own Discussion as prior
  context, then state clearly what is new: ff14SB/GBn2, off-lattice, over an enumerated space,
  with a certified optimum and a CI.**
* ~~**Still UNVERIFIED:** whether *Phys. Chem. Chem. Phys.* (2018), "Assessing AMBER force fields
  for protein folding in an implicit solvent" (PMID 29480910, doi 10.1039/c7cp08010g) covers
  **ff14SB + GB-Neck2**. RSC returned HTTP 403. Snippet-level only: *"the native state is not
  globally stable under the ff99SB force field in GB solvent."* **This is now the highest-value
  unfetched source in this document.**~~
  **⚠ RESOLVED, Sprint 16, 2026-09-06 (RETRACT): IT DOES.** Verbatim abstract obtained from the
  Europe PMC record. It tests **FF14SB and FF14SBonlysc with GB-Neck2**. **Credit where due: this
  bullet correctly identified the risk, correctly refused to assert it, and correctly pre-registered
  the consequence** ("Do it before any paper draft cites AMBER anti-ranking as novel"). Sprint 15 did
  **not** downgrade this source; it downgraded Maffucci & Contini and flagged this one. The reversal
  is of the *combined* §2.4 verdict, not of this bullet's judgement.

**Why this matters more than it looks.** If the peptide-folding literature already establishes
that ff14SB/GBn2 does not place the native at the global free-energy minimum for peptides, then
our "AMBER anti-ranks the native" result is **not a discovery about our pipeline — it is
documented expected behaviour of that force field in that regime.** That materially reframes
S14 novelty claim 3. **ACTION: this is the single highest-value primary fetch left undone in
this document.** Do it before any paper draft cites AMBER anti-ranking as novel.

---

## 3. HOW THE FIELD REPORTS RMSD — the section to cite when defending our protocol

### 3.1 RMSD is length-dependent, and the field knows it

> Carugo O, Pongor S. **"A normalized root-mean-square distance for comparing protein
> three-dimensional structures."** *Protein Science* **10**(7):1470-1473 (2001).

**[P15]** fetched. Verbatim: *"rmsd value obviously depends on the number of atoms included in
the structural alignment. Clearly, an rmsd value of, say, 3 A has a different significance for
proteins of 500 residues than for those of 50 residues."* They propose an `rmsd_100`
normalisation to a 100-residue reference.

**CAUTION — DO NOT QUOTE THE FORMULA.** The fetch returned **two mutually inconsistent
algebraic forms** for the normalisation within one response. **The exact `rmsd_100` formula is
UNVERIFIED.** Cite the qualitative length-dependence sentence (verbatim, uncontroversial) and
re-fetch before writing any algebra.

Roget et al. [P15] make the same point independently, more bluntly, and about our exact field:

> "An RMSD of 3 A is considered quite accurate for a protein of 100 amino acids but is awful
> for a peptide of 6 amino acids. On the other hand, quantum computing experiments tend to
> drastically downscale the problem (by using very small instances) and choosing specific
> instances used as proof of concept."

**This is the best single sentence in the literature for our Methods defence**, because it is
written by a quantum-PSP group about quantum-PSP practice, in 2026. It says both that our
3.204 A on 9-16mers is a *harder* number than 3 A on a 100-mer, and that downscaling to two
hand-picked miniproteins is a recognised failure mode of the field.

### 3.2 Six degrees of freedom in "the RMSD", all of them exercised in this literature

**Never place two RMSDs side by side without stating all six.**

1. **Atom set** — CA / backbone N,CA,C(,O) / heavy atoms / all atoms. QTF abstract says CA, and
   its Methods define nothing at all [S14]. Zhang et al. say "typically... backbone or Ca",
   undefined [P15]. Roget et al. say CA explicitly [P15]. **Ours: CA only.**
2. **Chain range** — full chain / terminal-trimmed / "structural core" / aligned region only.
   QTF states in a hardware aside that "the reported RMSD **excludes the terminal residues**
   and evaluates the structural core" [S14]. **On a 10-mer that is 20% of the chain**, and our
   own Sprint 12 measurement puts terminal dropout at **0.40-0.50 A cheaper than uniform
   dropout**. **Ours: all residues, no trimming.**
3. **Superposition** — Kabsch/optimal / none / TM-align / sequence-independent. **Ours: Kabsch,
   sequence-dependent.**
4. **Reference model** — crystal structure / one NMR model / ensemble medoid / lowest-energy NMR
   model. Gulsevin & Meiler use "the lowest-energy corresponding reference structure" [P15].
   **For 9-16mers most references are NMR ensembles and this choice is worth tenths of an
   Angstrom.**
5. **Normalisation** — raw / `rmsd_100` / length-normalised (the published Gulsevin & Meiler
   version reportedly uses "length normalized Ca RMSD" — UNVERIFIED). **Ours: raw.**
6. **Aggregation over the prediction set** — minimum over a pool (**ORACLE**) / pool mean /
   median of finals / the single selected model / best-of-5. **This axis carries the largest
   discrepancies in the literature and is the one most often left implicit.** **Ours: the single
   returned structure per target, then a paired mean across 126 targets with a bootstrap CI.**

### 3.3 The reporting standard we should adopt and defend

Emit, for every arm and every target set, the triple `s12/instrument.py` already produces:

    pool_best   (ORACLE DIAGNOSTIC — the representational / sampling ceiling)
    selected    (what the objective actually returns from the pool)
    returned    (what the full pipeline emits after aggregation and projection)

with the two gaps named. Roget et al. [P15] independently use the first two; QTF [S14] prints
both but leads its abstract with the oracle one; Gulsevin & Meiler [P15] effectively report the
third against the first. **No paper I found reports all three with the gaps named.** Doing so,
and arguing it should be mandatory in QPSP, is a cheap and defensible methodological
contribution. HYPOTHESIS that it is unclaimed; every component is LITERATURE-SUPPORTED.

---

## 4. THE COMPARABILITY MATRIX

**Our instrument, stated once, as the reference row for everything below:**

> **126 cluster-disjoint peptide targets, 9-16 residues. Full-chain Kabsch CA-RMSD over ALL
> residues, raw (not length-normalised), against the deposited reference. One returned
> structure per target. No oracle anywhere in the selection path. Paired mean across targets
> with a bootstrap 95% CI, five pinned leave-fold-out folds, a disjoint dev24, and a 60-target
> benchmark that is sealed until the architecture is frozen. Reproduces to 16 digits:
> `pool_best 1.7108244199364904 / top75_best 2.3061526409453816 /
> synthesis_fit 3.2040761603809194 / shipped 3.4540004952559396`.**

**Rule 1 (genuine VQE)**, as used in the `VQE?` column: a real parameterised quantum state, a
real variational optimisation of that state, and real sampling, with the objective evaluated as
a qubit observable. A classical optimiser searching real parameters through a
circuit-shaped nonlinearity onto a classically-scored structure is **NOT** genuine VQE by this
rule; it is a reparameterisation.

**Comparability verdict key:**
**[C]** directly comparable with stated caveats; **[I-oracle]** incomparable, the headline is an
argmin against the native; **[I-repr]** incomparable, a lattice/CG representation with a
different reachable set; **[I-metric]** incomparable, the RMSD definition differs or is
undefined; **[I-scale]** incomparable, n<=2 targets or a hand-picked instance; **[I-task]**
a different task entirely.

### 4A. Quantum protein / peptide structure prediction

| # | Paper | Problem | Representation | Dataset, lengths | Metric + exact definition | Info at prediction time | Predictive or ORACLE | Optimiser | Force field / energy | VQE genuine (rule 1)? | CVaR? | Samples | Major limitation | The exact gap | Comparability to our instrument | Ver |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **Q1** | **Cumbo et al., arXiv:2609.02113 (QTF), Sep 2026** | off-lattice all-atom PSP | continuous torsions from gauge-fixed statevector relative phases; all heavy atoms + side chains | **n=2**: chignolin 5AWL (10 res), Trp-cage 2JOF (20 res) | "Ca RMSD"; **no RMSD definition in Methods at all**; a hardware aside states it "excludes the terminal residues and evaluates the structural core" | sequence; **plus a hand-set E_e2e end-to-end bias whose target may be native-derived (UNVERIFIED)** | **BOTH, and the abstract leads with the ORACLE**: 0.623 / 2.501 A are argmins by RMSD; 2.90 / 5.39 A are energy-selected medians | COBYLA then SLSQP; best-of-50 scout init | custom hand-weighted; Rosetta REF15; **OpenMM AMBER ff14SB** | **NO.** Their words: the score "is not decomposed into quantum observables and measured as the expectation value of a qubit Hamiltonian" | **No.** Absent from the text and the 50-ref list | 7,334,433 structures; 1,200 replicas/protein; 8,192 shots/hardware job | no classical control; **zero CIs or p-values in 50 pages**; n=2; P ~ 2N+6n > N | no causal attribution of the quantum component | **[I-oracle][I-metric][I-scale]** on the headline. The 2.90/5.39 A medians are the only rows worth comparing, and even those are terminal-trimmed on 2 targets | [S14] |
| **Q2** | **Roget, Damour, Cadet, Wang, arXiv:2606.21241, Jun 2026** | **auditing the cost Hamiltonian, not predicting** | tetrahedral-lattice SAWs, CA only, arc 3.8 A | **12,446 peptides (4,866 unique seq), 6-600 aa**; complete enumeration at len<=15, 10,000 Pivot samples above | CA-RMSD after alignment vs experimental; explicitly defined | sequence only | **BOTH, and correctly separated**: min-RMSD labelled as the lattice ceiling; min-cost RMSD labelled as the achievable result | none (they enumerate/sample; no optimiser is run) | MJ contact energy, shell-resolved; HP-type as sensitivity check | N/A — no quantum run | No | full enumeration <=15; 10,000/peptide above (rho CI amplitude ~0.05) | lattice CA only; contact energies only | none for their own claim; **they identify the gap for everyone else** | **[I-repr]** for the numbers (~1.5 A lattice ceiling, ~5 A min-cost). **[C] for the METHODOLOGY** — this is our protocol, published | **[P15]** |
| **Q3** | **Zhang et al., arXiv:2510.06413, Oct 2025** | quantum + AI fused PSP | CA-only backbone from VQE; NetSurfP-3.0 SS/phi-psi/RSA priors | **75 fragments, 10-14 res, PDBbind**; 5 candidates each = 375 | "RMSD... typically computed over backbone or Ca atoms" — **atom set and chain range undefined** | sequence; NetSurfP-3.0 predictions | **PREDICTIVE.** `E_fuse` argmin; "The top-ranked structure c(1) is selected as the final prediction" | unnamed classical optimiser in the VQE loop | MJ pair energies + geometric/steric terms | Partially — VQE on IBM 127q with a Hamiltonian, but the fusion score is classical | **No** | 375 conformations; pool of **5** per target | pool of 5; PDBbind binding-site fragments; no fold/identity discipline; RMSD undefined | no held-out protocol, no CIs, tiny pool | **[C] with four caveats** — the closest predictive comparator that exists. **mean 4.89 A, median 4.70, sd 1.10** vs our 3.204 A | **[P15]** |
| **Q4** | **Kannan et al., arXiv:2510.15316, Oct 2025** | lattice PSP free-energy landscape | FCC lattice, turn encoding, 5(N-1) qubits | 13 peptides, 6-10 res (chains to 20) | RMSD vs experimental; definition not obtained | sequence | **ORACLE-flavoured**: the abstract reports **"minimum RMSD structures... 1.22 A to 3.11 A"** | CVaR-VQE + classical SA/MD controls | Miyazawa-Jernigan + non-bonded | Yes, on a lattice Hamiltonian | **YES** | not obtained | lattice; CG; MJ; the headline is a minimum | no predictive mean reported | **[I-oracle][I-repr]**. **This is the number our <2.0 A goal will be compared against, and it is a minimum over predicted structures on a lattice.** Not comparable | **[P15]** abstract, rest **[P13]** |
| **Q5** | **Yun, Seo, Jang, Park, Bae, Wu, arXiv:2606.01611, Jun 2026 (CD-QAOA)** | peptide PSP | tetrahedral lattice, CA | **n=1**: APRLRFY (7 res) | backbone RMSD, **against HF / DFT / MD structures, not against an experimental structure** | sequence | **PREDICTIVE** — "lowest energy from the probability distributions of quantum circuit outputs"; dominant bitstring `01001100110011`; **no native comparison at all** | COBYLA; then 10 ns explicit-water MD | growth-constraint + interaction Hamiltonian; MJ | Yes (CD-QAOA, lattice Hamiltonian) | No | l=1 and l=2; 12 qubits + 2 interaction qubits | **n=1 peptide, and the reference is a computed structure, not an experiment** | no experimental ground truth | **[I-scale][I-metric][I-task]**. RMSDs of 1.39/1.90/1.70 A are *to HF/DFT/MD*, not to a native. **Do not place beside any experimental RMSD** | **[P15]** |
| **Q6** | Robert, Barkoutsos, Woerner, Tavernelli, *npj QI* 7:38 (2021) | lattice PSP | tetrahedral turn encoding, 2 qubits/turn | 7-10-mers | native lattice fold recovery | sequence | predictive | CVaR-VQE, 20-qubit IBM | MJ-style contact | Yes | **YES — the canonical CVaR-QPSP paper** | — | lattice; CG; no all-atom | — | **[I-repr]**. Establishes CVaR-on-lattice-QPSP as taken | [P13] |
| **Q7** | Boulebnane, Lucas, Meyder, Adaszewski, Montanaro, *npj QI* **9**:70 (2023) | lattice peptide conformational sampling | lattice, self-avoiding walks + simplified potential | short peptides, 20-qubit simulations | energy/sampling quality, not RMSD | sequence | predictive | QAOA | simplified physical potential | Yes | No | 20 qubits, noiseless | simplified setting by their own account | — | **[I-task]** for numbers; **the mandatory control paper** for methodology | **[P15]** abstract |
| **Q8** | Pamidimukkala et al., *JCTC* 20:10223 (2024) | gate-based PSP, high DOF | turn-based with greater DOF; HP alphabet | varied lengths, up to 114 qubits on IBM | folded structures | sequence | not established | VQE | HP | Yes | not established | — | HP simplification; linear, large qubit count | — | **[I-repr]** | [S] |
| **Q9** | Li, Doga, Raubenolt, ... Blankenberg, arXiv:2507.08955 (Jul 2025) | FCC-lattice PSP with constraints | FCC lattice | **n=1**: KLVFFA (6-mer) | — | sequence | — | VQE-C, Lagrangian duality; IBM Eagle R3 + Heron R2 | contact potential | Yes | No | — | single 6-mer; lattice | — | **[I-scale][I-repr]** | [S14] abstract |
| **Q10** | Linn, Li, ... Blankenberg et al., arXiv:2509.18263 (Sep 2025) | problem-agnostic-ansatz PSP | tetrahedral / BCC / FCC, 2nd-nearest-neighbour | peptides to 26 aa | — | sequence | — | hardware-efficient problem-agnostic ansatz; IBM Kingston | "classically computable energy cost function" | **Partially** — the cost is computed classically, same house style as QTF | No | — | lattice; results not in abstract | — | **[I-repr]** | [S14] abstract |
| **Q11** | Zhang et al., QSAD, arXiv:2607.06971 (Jul 2026) | binding-pocket peptide reconstruction | amino-acid-level Hamiltonian sampling, non-iterative | **101 peptides, 5-18 res** | RMSD after optimal superposition | sequence + pocket context | not established | non-iterative Hamiltonian evolution; IBM Heron R2 | — | Partially | No | — | binding-pocket regime; improvement stated only in relative terms ("27-71%") | absolute accuracy not extractable | **[I-task]**, but **the only quantum peptide benchmark at our scale.** Watch it | [S14] abstract |
| **Q12** | arXiv:2506.22677 / *Adv Sci* (2026) | binding-site structure prediction | tetrahedral lattice, sparse Pauli Hamiltonian with steric/geometric/chirality terms | binding sites | — | sequence | — | VQE, two-stage (energy estimation separated from measurement decoding) | problem-specific Hamiltonian | Yes | No | — | binding-site niche | — | **[I-task]** | **[P15]** snippet |
| **Q13** | Feng et al., Graph-VQE, arXiv:2607.02749 (Jul 2026) | partitioned-Hamiltonian PSP | lattice-style interaction graph | not specified | **reports ENERGY, not RMSD** | — | — | Louvain graph partitioning; simulated multi-QPU CUDA-Q | — | Yes | No | — | **never prices the structure** | — | **[I-task]**. A live example of the failure our program names | [S14] abstract |
| **Q14** | Casares, Campos, Martin-Delgado, QFold, arXiv:2101.10279 / *QST* 7:025013 | off-lattice torsion PSP | **torsion angles, user-defined discretisation over the full 0-360 range** | small peptides | — | sequence + Minifold dihedral initialiser | — | quantum-walk Metropolis; minimal IBMQ Casablanca realisation | Psi4 energies | Different paradigm (quantum walk, not VQE) | No | — | tiny systems; angles binned to few bits | — | **[I-repr]**. **Kills the unqualified "first off-lattice quantum PSP" claim** | **[P15]** search-level; primary still not fetched |
| **Q15** | Perdomo-Ortiz et al., *Sci Rep* 2:571 (2012) | lattice HP folding | 2D/3D lattice | 6-mer | lattice ground state | sequence | predictive | D-Wave annealing | HP | N/A | No | — | toy lattice | — | **[I-repr]** | [S] |
| **Q16** | Fingerhuth, Babej, Ing, arXiv:1810.13411 (2018) | QAOA lattice folding | lattice turns | lattice peptides | — | sequence | — | quantum alternating operator ansatz | HP/MJ | Yes | No | — | lattice only | — | **[I-repr]** | [S] |
| **Q17** | Chandarana et al., *Phys Rev Appl* 20:014024 (2023) | digitised counterdiabatic QAOA | lattice turns | lattice peptides | — | sequence | — | DC-QAOA | HP/MJ | Yes | No | — | lattice | — | **[I-repr]** | [S] |
| **Q18** | Marchand et al. (2018) | torsion-space quantum optimisation | rotatable-bond torsions | small molecules | — | — | — | D-Wave annealing | **"flexible with respect to the choice of molecular force field"** | N/A | No | — | small molecules, not peptides | — | **[I-task]**. **Torsion-space quantum optimisation with a real force field, 8 years before QTF** | [P13] |
| **Q19** | Mato et al., Quantum Molecular Unfolding (2022) | torsion HUBO | torsions over rotatable bonds | ligands | — | — | — | annealing / QAOA | force field | N/A | No | — | ligand scale | — | **[I-task]** | [P13] |
| **Q20** | Agathangelou, Manawadu, Tavernelli, arXiv:2507.19383 (2025) | side-chain packing | fixed backbone, rotamer choice | side-chain packing | — | backbone given | — | VQE / QAOA **with an explicit classical comparison** | rotamer energy | Yes | not established | — | fixed backbone | — | **[I-task]**, but **the only quantum-vs-classical comparison in this table.** A model for how to run ours | [S] |
| **Q21** | Doga, Raubenolt, Cumbo et al., *JCTC* 20(9):3359 (2024); arXiv:2312.00875 | perspective / problem-selection framework | survey | Zika NS3 helicase catalytic loop as proof of concept | — | — | — | — | — | — | — | — | a perspective, not a benchmark | — | **[I-task]**. The field's own framing document; same group as QTF | **[P15]** search-level |
| **Q22** | Shajan, Kaliakin, Merz et al. (SQD / SQD+DMET), 2025-26 | all-atom quantum **chemistry** on proteins | electronic structure | protein-ligand complexes to 12,000 atoms | electronic energies | — | — | SQD | — | genuine quantum chemistry | No | — | **scores conformers; does not predict folds** | — | **[I-task]**. The other meaning of "quantum protein" | [S] |

### 4B. Variational quantum algorithms — objective, geometry, trainability, controls

| # | Paper | Contribution | Relation to our claims | Ver |
|---|---|---|---|---|
| V1 | **Barkoutsos, Nannicini, Robert, Tavernelli, Woerner. *Quantum* **4**:256 (2020). doi:10.22331/q-2020-04-20-256. arXiv:1907.04769** | **CVaR as the VQE aggregation function** for diagonal (classical) Hamiltonians; "CVaR as an aggregation function leads to faster convergence to better solutions for all combinatorial optimization problems tested" | The origin of our CVaR arm. **The paper contains no random-sampling or untrained-circuit baseline** (checked at the journal page level; the fetch returned no such comparison, **UNVERIFIED as an absence**) | **[P15]** journal page; full PDF not read this sprint |
| V2 | Kandala et al., *Nature* 549:242 (2017) | hardware-efficient VQE; the EfficientSU2 lineage | QTF's ansatz; our `StatevectorCircuit` family | [S] |
| V3 | McClean, Boixo, Smelyanskiy, Babbush, Neven, *Nat Commun* 9:4812 (2018) | barren plateaus from 2-design expressibility | Scope: needs the circuit ensemble near a 2-design. Our MPS/shallow arms are outside it | [P13] |
| V4 | Cerezo, Sone, Volkoff, Cincio, Coles, *Nat Commun* 12:1791 (2021) | cost-function-dependent BPs; global observables decay even at shallow depth | Our measurement refutes its relevance at 6-18 qubits: the measured ansatz kernel is flat in Pauli weight | [P13] |
| V5 | Anschuetz & Kiani, *Nat Commun* (2022) | **traps, not plateaus**: BP-free shallow models still have superpolynomially few good minima | **Absence of a plateau is not evidence of trainability.** Directly supports making solution quality its own axis | [P13] |
| V6 | Qiu, Lumbreras, Li, Rebentrost, arXiv:2605.02850 (2026) | the tilted-loss family (containing CVaR) "does not remove the barren plateau problem by itself"; the bottleneck moves "from trainability to estimability" | **CVaR must never be justified as BP mitigation.** Justify it as a tail objective | [P13] |
| V7 | Cerezo et al., *Nat Commun* (2025), simulability critique | many provably BP-free models admit classical simulation after a classical data-acquisition phase | Live against our exact-statevector and exact-MPS arms. Unaddressed by us and by QTF | [P13] |
| V8 | Ragone et al. (2024), Lie-algebraic variance unification | `Var = P_M(rho)*P_M(O)/dim(M)` | The language for our spectrum-to-gradient-variance chain | [P13] |
| V9 | **Kang P., arXiv:2605.01319 (May 2026, v2 Aug 2026), "Barren Plateaus as Destructive Interference" / "Decomposing Gradient Suppression in Barren Plateaus"** | **exact term-resolved decomposition of the gradient second moment into pre-cancellation activity, sign organisation, and coupling** | **The closest published relative of our spectrum-to-variance chain.** Kang decomposes by Hamiltonian TERM; we predict the variance from the Walsh/Pauli SPECTRUM of a molecular energy with no free parameters. Adjacent, not the same. **Fetch before finalising claim 2** | **[P15]** search-level |
| V10 | Fourier-expansion-of-VQA literature (e.g. *Phys Rev A* 108:032406) | the VQA loss is a partial Fourier series; coefficients computable classically to degree m in O(N 2^m); coefficient variance proportional to frequency redundancy | Establishes that spectrum-to-variance reasoning is a known *technique*. **Our claim must be scoped to applying it to a molecular energy, not to inventing it** | **[P15]** search-level |
| V11 | Stokes, Izaac, Killoran, Carleo — quantum natural gradient; Fubini-Study / QFIM geometry | QNG as steepest descent under the Fubini-Study metric | Our "QNG is refuted ONLY at depth 1; open again at depth >= 2 (condition number 155.2 at d=2, 478.9 at d=3)" sits inside this literature. **Whether QNG HELPS at depth >= 2 for a molecular objective is our open HYPOTHESIS** | **[P15]** search-level |
| V12 | Mariyanto et al., *Int J Quantum Chem* (2026), doi 10.1002/qua.70255 | Fubini-Study metric tensor evolution and reduced-state distinguishability in geometry-aware VQE optimisation | Directly relevant to the depth>=2 QNG question. **Fetch before running a QNG arm** | **[P15]** search-level |
| V13 | VQEC — Purdue/Kekatos et al., arXiv:2311.08502; *Phys Rev A* 110:022430 (2024) | constrained VQE via a Lagrangian over circuit parameters AND dual variables, perturbed primal-dual with the parameter-shift rule | The mature form of constrained VQE. **Our binary encoding is surjective and needs no penalty**, which sidesteps this whole line — a genuine simplification worth stating | **[P15]** search-level |
| V14 | Augmented-Lagrangian / slack-free / penalty-free quantum optimisation (arXiv:2503.10077, 2507.12159, 2301.12393) | penalty parameters need not diverge; slack-free reformulations save qubits; penalty-based optimisers "problematic and sometimes unable to properly converge" | Corroborates our S14 encoding finding: **a non-surjective encoding's Pauli spectrum measures its CONSTRAINT, not its objective** | **[P15]** search-level |
| V15 | Quantum circuit Born machines / MPS Born machines / tensor-network + autoregressive quantum states | generative quantum models; MPS-based Born machines support tractable log-likelihood and autoregressive sampling | **The natural alternative generator to a VQE ansatz.** Relevant to our `MPSAnsatz`. No QPSP application found | **[P15]** search-level |
| V16 | Hybrid CV-DV variational quantum computing (arXiv:2608.03907, 2606.05297, 2511.13882, HyQBench arXiv:2603.04398) | bosonic modes represent continuous variables directly, "avoiding the qubit overhead associated with discretization and binary encoding" | **The principled way to do continuous torsions on quantum hardware** — and the correct contrast to QTF's phase encoding. Not applied to PSP anywhere I found | **[P15]** search-level |
| V17 | Riemannian / Stiefel-manifold quantum circuit optimisation (arXiv:2202.06976, 2501.08872, 2506.17395, 2602.20605) | gradient flow on SU(2^n) / the complex Stiefel manifold; exact geodesic transport | An alternative to QNG for the depth>=2 conditioning problem. Untested on a molecular objective | **[P15]** search-level |
| V18 | Truchon & Bayly, *JCIM* 47:488 (2007) — BEDROC / RIE; and CROC (*Bioinformatics* 26:1348, 2010) | **the "early recognition" problem**: ROC/AUC/average-rank are the wrong metrics when only the top of a ranked list matters; BEDROC with alpha=20 puts 80% of the weight on the top 8% | **This is the prior art for tail-restricted discrimination analysis.** The *concept* that top-of-list skill differs from bulk skill is standard in virtual screening. **Our claim must be scoped to the specific finding (objectives are at chance INSIDE their own tail, against a random-tail null at 0.524-0.527), not to the idea of tail-restricted evaluation** | **[P15]** search-level |

### 4C. Classical structure prediction, scoring, force fields, geometry

| # | Paper | What it establishes | Metric definition, where it matters | Relation to us | Ver |
|---|---|---|---|---|---|
| C1 | **Handl, Knowles, Lovell, *Bioinformatics* 25(10):1271-1279 (2009)** | decoy-set evaluation is artefact-ridden; native is an extremum of an Amber99 term in **139/149** sets; rank/z-score of the native is a weak test | native-rank / z-score | **Confirms our "decoy-bank ranking does not transfer". Our finding is a replication, not a discovery** | **[P15]** |
| C2 | **Gulsevin & Meiler, bioRxiv 2022.02.17.480937 / *Structure* 31(1) (2023)** | AF2 on peptides; per-class mean CA-RMSD 2.2-5.7 A; **only 13% of lowest-RMSD models were AF2's own rank-0**, vs a 20% uniform null | CA-RMSD, CA alignment, vs the lowest-energy NMR reference model; preprint says 155 peptides 16-60 aa (journal version reportedly 588, 10-40 aa, length-normalised — **UNVERIFIED**) | **The state of the art fails our exact selection test, on peptides.** The single best "nothing ranks within the pool" citation | **[P15]** |
| C3 | **Carugo & Pongor, *Protein Science* 10(7):1470-1473 (2001)** | RMSD is length-dependent; proposes `rmsd_100` | **the formula returned inconsistently on fetch — UNVERIFIED, do not quote it** | The citation for why a 3 A on a 12-mer is not a 3 A on a 100-mer | **[P15]** |
| C4 | **Maffucci I, Contini A. *JCTC* 12(2):714-727 (2016), doi:10.1021/acs.jctc.5b01211** | REMD over **six AMBER force fields x three implicit solvent models** on 2 helical + 3 beta-hairpin + 3 IDP peptides. Verbatim: **"a combination of the force field and implicit solvation models able to accurately predict the native structure of all the considered peptides was not identified."** **Force fields tested are ff96 and the ff99SB/ildn/ildn-phi series — NOT ff14SB** | secondary-structure agreement, not RMSD ranking | Prior context for AMBER-family peptide failure, but **it does not document our ff14SB/GBn2 pair**, so the threat to novelty claim 3 is **downgraded**. Cite it in the Discussion, then state what is new | **[P15]** abstract |
| C4b | **Shao Q, Zhu W. *PCCP* 20(10):7206-7216 (2018), doi:10.1039/c7cp08010g, PMID 29480910** | **RESOLVED, Sprint 16 2026-09-06 (RETRACT): IT DOES COVER OUR PAIR.** Enhanced-sampling MD over six AMBER force fields including **FF14SB and FF14SBonlysc**, each with **GB-Neck2**, on TC10b (20 res, α), HP35 (35, α), **1E0Q (17 res, β-hairpin)** and GTT (~35, β). Verbatim: *"the measured counterparts are significantly discrepant in the cases of larger or β-structured peptides (HP35, 1E0Q, and GTT)"*; *"a combination ... able to describe all aspects of the folding transitions towards the native structures of all the considered peptides was not identified"*; but also *"FF14SBonlysc coupled with the GB-Neck2 model seems to be a reasonably balanced combination"* | **folding thermodynamics / free-energy profiles from enhanced-sampling MD** — NOT conformer ranking of a fixed pool, and no rank statistic appears | ~~downgraded~~ **REINSTATED [S]-tier threat, RE-SCOPED.** "AMBER anti-ranks peptide natives" is documented expected behaviour, not a discovery. Claim only the **observable** (rank position vs a certified optimum), the **statistics**, the **AMOEBA extension**, and the **stereochemical-repair** result. The length defence FAILS (1E0Q is 17 res); the "free energy is weaker" defence is BACKWARDS (it is stronger). §2.4 correction box | **[A+]** publisher abstract record via Europe PMC; full text still not obtained |
| C5 | Maier, Martinez, Kasavajhala, Wickstrom, Hauser, Simmerling, *JCTC* 11:3696 (2015) | AMBER **ff14SB** | — | The identical force field we run via `core/amber.py`; also QTF's OpenMM backend | [S] |
| C6 | Eastman et al., *PLoS Comput Biol* 13:e1005659 (2017) | OpenMM 7 | — | The identical engine we run | [S] |
| C7 | Alford et al., *JCTC* 13:3031 (2017) | Rosetta REF15 | — | QTF's second backend; **its final-model native-like hit rate on chignolin was 0.3%** | [S] |
| C8 | Lindorff-Larsen et al., *Proteins* 78:1950 (2010) | amber99sb-ildn | — | QTF's GROMACS post-minimisation force field | [S] |
| C9 | Miyazawa & Jernigan, *J Mol Biol* 256:623 (1996); *Macromolecules* 18:534 (1985) | MJ contact potentials | — | The energy of nearly every row in 4A. **Roget et al. [P15] note MJ was parameterised on proteins >50 residues — which is why it anti-correlates on our length regime** | [S] |
| C10 | Dill (1985); Lau & Dill (1989) | HP model | — | The other lattice energy | [S] |
| C11 | Zhou & Zhou, DFIRE, *Protein Sci* 11:2714 (2002) | distance-scaled statistical potential | — | Closest published analogue of our Legacy energy | [S] |
| C12 | Parsons, Holmes, Rojas, Tsai, Strauss, *J Comput Chem* 26:1063 (2005) — NERF; AlQuraishi pNERF, *J Comput Chem* 40:885 (2019) | torsion -> Cartesian reconstruction | — | QTF's reconstructor; ours is `core/geometry.build_backbone_batch` on ideal geometry — same class | [S] |
| C13 | Engh & Huber, *Acta Cryst* A47:392 (1991) | ideal bond lengths and angles | — | **Shared assumption, shared ceiling** with QTF | [S] |
| C14 | Honda et al., *JACS* 130:15327; PDB 5AWL. Barua et al., *PEDS* 21:171; PDB 2JOF | chignolin; Trp-cage | — | QTF's two targets. **Both are among the most-optimised miniproteins in structural biology** | [S] |
| C15 | Coarse-grained protein model literature (*Chem Rev* 116:7898, 2016; *PLoS Comput Biol* 6:e1000827) | CG landscapes: "non-native structures... with energies similar to those of the native state"; difficulty "for small peptides" | — | The classical statement of our Legacy problem | [S] |
| C16 | CASP model-quality / consensus QA literature (QMEANclust, *BMC Struct Biol* 9:35, 2009, and successors) | consensus QA "tend[s] to fail if the best models are far from the dominant structural cluster" | — | **The external form of "consensus is outlier avoidance, capped at the pool mode"** | **[S]** |
| C17 | ML scoring-function generalisation: *JCIM* 2023 doi 10.1021/acs.jcim.2c01149; arXiv:2512.05386; arXiv:2605.11764 | random-CV > sequence-CV > Pfam-CV monotone degradation; "random splits reward within-target interpolation, whereas leave-one-target-out measures the novel-target prediction"; the gap is **"documented but not decomposed"** (arXiv:2605.11764) | — | **The within-target/cross-target collapse is DOCUMENTED for protein-LIGAND affinity.** Our claim must be scoped to a *structural conformational* objective. See §5 | **[P15]** search-level |
| C18 | Homology-leakage / benchmark-contamination literature (FoldBench, *Nat Commun* 2025; *Brief Bioinform* 26(2):bbaf104, 2025; SafeBench-Seq arXiv:2512.17527) | cluster-level holdout at <=40% identity is the accepted standard; benchmark performance overestimates generalisation without it | — | **Directly supports our pinned identity clusters and leave-fold-out discipline.** Our memory entry "benchmark and folds must be pinned" is the field standard, not an idiosyncrasy | **[S]** |
| C19 | Cornilescu, Delaglio, Bax, TALOS, *J Biomol NMR* 13:289 (1999); Shen et al. TALOS+, 44:213 (2009); **Shen & Bax, TALOS-N, *J Biomol NMR* 56:227 (2013)** | shifts -> (phi,psi); TALOS-N: **>=90% of residues, error rate < ~3.5%, RMS ~12 deg vs crystallographic (phi,psi), chi1 rotamer for about half** | — | **Re-confirmed at abstract level this sprint [P15].** Channel already closed by arithmetic in S14 (ORACLE-perfect torsions on all 54 covered targets leave the instrument at 2.021 A) | **[P15]** abstract + [P13] |
| C20 | Shen & Bax, SPARTA+, *J Biomol NMR* 48:13 (2010); Li & Bruschweiler, UCBShift (2020); LEGOLAS, *JCTC* (2025) doi 10.1021/acs.jctc.5c00026 | coordinates -> shifts (the inverse map) | — | Inverse maps, not torsion predictors. **No transformer-based shift->torsion successor to TALOS-N found in a targeted 2025-2026 search.** LITERATURE-SUPPORTED as an absence | [P13] + **[P15]** absence |
| C21 | Klukowski, Riek, Guntert, ARTINA, *Nat Commun* 13:6151 (2022); *Sci Adv* 9:adi9323 (2023) | raw NMR spectra -> structure; 1.44 A median RMSD, 91.4% assignment on a 100-protein benchmark | — | **Numbers still snippet-level. Re-fetch before pricing the channel** | [S] |
| C22 | SPOT-1D-LM, *Sci Rep* (2022) | sequence -> (phi,psi); phi MAE 15.99 deg / psi 23.74 deg on TEST2018 | — | **No 9-16-mers in any test set**; a 13-mer is the Neff~1 regime | [P13] |
| C23 | Distance-geometry / NMR restraint literature (X-PLOR-NIH DG; metric-matrix embedding) | DG "is sensitive to small uncertainties in the distance matrix"; extended polypeptides are "usually underdetermined... neither distance geometry nor ab initio simulated annealing will produce unique structures" | — | **Directly supports our memory entry "the distance prior is the ceiling."** Underdetermination of short extended chains is a known DG result | **[S]** |
| C24 | Peptide conformational-ensemble literature (accelerated MD / REMD benchmarking, e.g. PMC10207344) | short peptides are conformationally heterogeneous; single-structure comparison is questionable for them | — | Supports treating our 9-16mers as ensembles and reporting the sampled-vs-selected gap | **[S]** |

**Row count: 22 (4A) + 18 (4B) + 24 (4C) = 64 rows; the BRIEF asks for >= 40.**
**Primary-verified this sprint [P15]: Q2, Q3, Q4(abs), Q5, Q7(abs), Q12(snip), Q14(snip), Q21(snip), V1, V9-V18(search), C1, C2, C3, C17, C19, C20.**

### 4D. The comparability verdict, in one paragraph

Of 22 quantum rows, **one (Q3) is even indicatively comparable to our instrument**, and it
carries four caveats. Two more (Q1's medians, Q2's methodology) are usable with explicit
adjustment. **Every remaining headline accuracy number in quantum PSP is incomparable to a
126-target, full-chain, no-oracle CA-RMSD mean**, for at least one of four reasons: it is an
argmin against the native (Q1, Q4, and by construction most lattice work), it lives on a
lattice whose reachable set differs (Q4, Q6, Q8-Q10, Q15-Q17), its RMSD is undefined or
terminal-trimmed (Q1, Q3), or it is one or two hand-picked targets (Q1, Q5, Q9). **This is
itself the strongest argument for our protected instrument, and it should be a table in the
paper.**

---

### 4E. THE CVaR-VQE PEPTIDE SUB-LITERATURE THAT SPRINT 14 MISSED

This is the second claim-threatening discovery of the sprint and it deserves its own block.
S14 concluded "genuine in-loop CVaR remains unclaimed in QPSP torsion space" on the basis that
QTF has no CVaR and that Robert (2021) / Kannan (2025) are lattice work. **That is still true
of torsion space, but the lattice side is far more developed than S14 recorded.** There is a
standing group publishing CVaR-VQE peptide folding repeatedly.

| # | Paper | What it is | Ver |
|---|---|---|---|
| **X1** | **Uttarkar A, Niranjan V, Saxena A, Kumar V. "QuPepFold: A python package for hybrid quantum-classical protein folding simulations with CVaR-optimized VQE." *PLoS One* 21(2):e0342012 (2026). doi:10.1371/journal.pone.0342012** | **Tetrahedral lattice**, turn encoding, `2(N-3)` turn qubits + interaction qubits. **CVaR in the optimisation loop**, `CVaR_alpha(theta) = (1/alpha) sum_{i in S_alpha} p_i E_i`, **alpha = 0.025 noiseless / 0.05 on hardware**. Energy `H = lambda_dis*E_dist + lambda_loc*E_local + lambda_back*E_backfold + H_MJ`, modified Miyazawa-Jernigan, no solvation. **1,224 sequences, 6-10 aa, ~1.76 million conformers.** CVaR-VQE reaches the ground state **~30% faster** than expectation-value VQE. **Reports NO RMSD at all** — only energies (kcal/mol), success rates, and unique-geometry counts. No classical folding control | **[P15]** full text fetched |
| **X2** | Uttarkar et al. (same group). "Quantum synergy in peptide folding: A comparative study of CVaR-variational quantum eigensolver and molecular dynamics simulation." *Int J Biol Macromol* (2024), PMID 38862055 | **50 peptides of 7 residues.** CVaR as the aggregation function, **100 iterations x 500,000 shots**, against 50 ms MD. Concludes CVaR-VQE gives "more effective folding outcomes with respect to sampling and global optimization" | **[S]** snippet only |
| **X3** | Same group. "A comparative insight into peptide folding with quantum CVaR-VQE algorithm, MD simulations and structural alphabet analysis." *Quantum Information Processing* (2024), doi:10.1007/s11128-024-04261-9 | CVaR-VQE + MD + structural alphabet analysis | **[S]** snippet only |

**Three consequences.**

1. **"In-loop CVaR for peptide folding" is TAKEN.** Not narrowly, not once — by a group with at
   least three papers, a released package, and an explicit in-loop `CVaR_alpha` with stated
   alpha values that bracket ours. S14's phrasing ("CVaR-VQE is taken on lattices") was right
   but understated the depth.
2. **What survives is much narrower than S14 stated**: CVaR in a **continuous / high-resolution
   off-lattice torsion representation with a real all-atom force field**, and — more durably —
   **the CVaR defect analysis itself.** See §5 claim 4.
3. **X1 is a fourth instance of the field-wide pathology.** A 2026 *PLoS One* paper on peptide
   folding, 1.76 million conformers, reports **no RMSD to any experimental structure**. It
   optimises the objective and never prices the structure — exactly Graph-VQE (Q13) again.
   **This strengthens the negative-evidence story considerably: it is not that the field gets
   bad RMSDs, it is that a substantial fraction of it does not compute them.**

---

## 5. THE NOVELTY BOUNDARY — each S14 claim re-tested against the current literature

Verdicts: **SURVIVES** / **SURVIVES (RE-SCOPED)** / **TAKEN** (with citation) / **UNVERIFIED**.
I have tried to break each one. Where I broke it, I say so.

### 5.1 The six claims S14 said survive

#### Claim 1 — the exact torsion-space locality theorem. **SURVIVES.**

`d_ij` depends on exactly the `j-i-1` residues strictly between i and j (agreement 1.0000);
all-atom supports one residue wider (CA `i<m<j`, N `i<=m<j`, C/O `i<m<=j`, CB `i<=m<=j`);
corollary that a separation-8 pair is a 14-qubit interaction at k=4, so **no 2-local Ising form
of a distance-based molecular objective exists in this encoding**.

**Re-tested.** A targeted search for torsion-space locality / exact support of pairwise
distances in internal coordinates returned only (i) the standard observation that fixed bond
geometry leaves torsions as the DOF, and (ii) pairwise *statistical* dependence of torsions in
cheminformatics — a different object entirely. **Nothing states or proves the support result.**
Roget et al. [P15] do not analyse locality; QTF [S14] does not, and does not need to because it
never forms a qubit Hamiltonian. **SURVIVES, and it remains our strongest formal claim.**
Caveat: a negative literature search is weaker evidence than a positive one; the result is
elementary enough that it may exist in the robotics / molecular-kinematics literature under
different language (kinematic chains, Denavit-Hartenberg). **HYPOTHESIS: worth one search pass
in the robotics literature before publication.** That search is NOT done.

#### Claim 2 — the no-free-parameter spectrum-to-gradient-variance chain. **SURVIVES (RE-SCOPED).**

Measured/predicted 1.006 and 1.001, plus the delta-spike artefact correction (raw AMBER's
top-10 of 4,096 configurations carry a median 99.6% of Walsh variance, landing on
`Binomial(m,1/2)` exactly; 99th-percentile winsorisation insufficient).

**Re-tested, and the scope must shrink.** Two adjacent lines exist:
* **Fourier/Walsh analysis of VQA loss functions is an established technique** (V10): the loss
  is a partial Fourier series, coefficients are classically computable to degree m in
  `O(N 2^m)`, and coefficient variance is proportional to frequency redundancy.
* **Kang, arXiv:2605.01319 (2026)** (V9) gives an *exact* decomposition of the gradient second
  moment into activity, sign organisation and coupling — a term-resolved variance theory.

**Therefore: "spectrum predicts gradient variance" is NOT novel as a principle.** What appears
unclaimed is the specific instantiation: **computing the Walsh/Pauli spectrum of a real
molecular energy (AMBER ff14SB and a coarse energy) over an enumerated configuration space, and
predicting the measured gradient variance from it with no fitted parameter.** Nothing in the
QPSP literature computes a Pauli spectrum of a molecular energy at all. **SURVIVES as an
application and as a measurement; does NOT survive as a theoretical device.** The delta-spike
artefact and the winsorisation-insufficiency correction appear genuinely new and are the more
defensible half. **ACTION: fetch Kang arXiv:2605.01319 in full before drafting this claim.**

#### Claim 3 — the certified global optimum on an enumerable space. **TAKEN (RE-SCOPED to survive).**

> **TAKEN BY: Roget M, Damour C, Cadet F, Wang J. arXiv:2606.21241 (19 June 2026).**

They enumerate the complete feasible space for every peptide of length <= 15 across 12,446 PDB
peptides, compute the RMSD of the cost function's certified global minimum, and state that
**"the minimum-cost conformation has, on average, a larger RMSD than a randomly selected
feasible conformation."** That is our result, in print, three months ago.

**What survives, precisely:**

* Their space is a **tetrahedral CA lattice**; ours is an **off-lattice sequence-conditioned
  torsion library at k=4** with real backbone geometry. Different reachable sets, different
  claim.
* Their objective is a **Miyazawa-Jernigan contact potential** (plus an HP sensitivity check);
  ours are **AMBER ff14SB/GBn2 all-atom via OpenMM** and a non-all-atom Legacy energy. **No one
  has published a certified global optimum of an all-atom force field over an enumerated
  conformational space.** That is because nobody else has both an enumerable space and 28 ms
  all-atom single points.
* Their measurement is a **sign statement** ("worse than random on average"). Ours is a
  **magnitude with a paired CI** (+0.139 A, against a matched random arm), plus the
  decomposition of *why* (filtering vs ordering), plus the tail analysis.

**Verdict: SURVIVES ONLY AS "the certified global optimum of an all-atom force field, off
lattice, quantified with a CI and decomposed." The bare phenomenon is TAKEN and must be cited
to Roget et al., not claimed.** This is the single largest change to the novelty story since
S14, and the coordinator should treat it as binding.

**Bonus, and it cuts the other way — they hand us a mechanism.** Their length curve (rho
negative for short peptides, ~0 at length 50, positive above) with the explanation that **MJ
was parameterised on proteins longer than 50 residues** suggests our energy-ranking failure may
be a **length-regime artefact of potentials fitted on globular proteins**, not a universal
property. **HYPOTHESIS, testable, and it is the most interesting new research question in this
document.** Our AMBER arm is *not* a knowledge-based potential fitted on large proteins, so if
AMBER also anti-ranks at 9-16 residues, their explanation does not cover us — and *that*
contrast is publishable.

#### Claim 4 — in-loop CVaR in torsion space. **SURVIVES, VERY NARROWLY, AND THE GROUND MOVED.**

**Newly TAKEN this sprint:** in-loop CVaR-VQE for peptide folding, on a lattice, with stated
alpha values — **Uttarkar, Niranjan, Saxena, Kumar, *PLoS One* 21(2):e0342012 (2026)**
(QuPepFold: `alpha = 0.025` noiseless, `0.05` hardware, in-loop, 1,224 sequences 6-10 aa), plus
two earlier papers from the same group (X2, X3), on top of Robert et al. (2021) and Kannan et
al. (2025).

**What still survives:**
* **CVaR in a continuous / high-resolution off-lattice torsion representation, scored by an
  all-atom force field.** Not found anywhere. QTF has no CVaR at all.
* **Far more durably: the CVaR defect analysis itself.** Nothing in the CVaR literature I
  searched reports (i) the closed form of the `baseline="tail"` gradient bias
  `-c * grad P(E<q)`, (ii) the sampled-CVaR upward bias at non-integer `alpha*N` (+0.134 sd at
  N=13, alpha=0.1, decaying as 1/N), or (iii) **`dCVaR/dp` identically zero iff
  `p(argmin E) >= alpha`, with 0/3000 counterexamples, so that at alpha=0.01 one initialisation
  in four to twenty starts dead.** Barkoutsos et al. [V1] introduce CVaR and report faster
  convergence; they do not analyse its estimator or gradient pathologies. **These three defects
  are the most defensible unclaimed CVaR result we hold, and they are stronger than "we used
  CVaR in torsion space."** Reposition accordingly.
* The correction that **small alpha discriminates better and concentration collapses variants**
  (the recorded claim had the sign backwards) is a methodological correction of independent value.

#### Claim 5 — the causal VQE-vs-classical control at matched budget. **SURVIVES.**

**Re-tested hard.** The comparison landscape:
* **Boulebnane et al. (2023)** [Q7] compare QAOA against **uniform random sampling** and find
  parity "up to a small overhead". Closest precedent. **Not the same control.**
* **Agathangelou/Tavernelli (2025)** [Q20] compare quantum against classical for **side-chain
  packing on a fixed backbone.** Different problem.
* **Kannan et al. (2025)** [Q4] compare against classical simulated annealing and MD, but
  report a minimum RMSD, not a matched-budget paired comparison.
* **Uttarkar et al. (X1)** compare CVaR-VQE against **expectation-value VQE** — an internal
  ablation, not a classical control — and explicitly provide "no head-to-head benchmarking of
  their package against classical folding tools."
* **QTF** [Q1] has none, over 50 pages.

**SURVIVES.** A matched-evaluation-budget, paired, CI-reported VQE-vs-classical control on a
molecular objective with a defined RMSD protocol is not in this literature. Per BRIEF §1.6, a
rigorous negative here is a legitimate success.

#### Claim 6 — the protected statistical instrument. **SURVIVES, and it is the cheapest to defend.**

QTF: **zero CIs or p-values in 50 pages, n=2 targets** [S14]. Q5 (CD-QAOA): n=1 peptide, and its
reference structures are **HF/DFT/MD models, not experiments** [P15]. Q3: no fold discipline,
pool of 5 [P15]. Q4: a minimum over 13 peptides [P15]. X1: 1,224 sequences and **no RMSD at
all** [P15]. The only row with real statistical scale is Q2 (Roget et al., 12,446 peptides) —
and **Q2 runs no predictor**; it is an audit.

**Corroborated externally**: cluster-level holdout at <=40% identity is the accepted
anti-leakage standard (C18), so our pinned identity clusters and leave-fold-out folds are field
standard rather than idiosyncratic. **SURVIVES.** Nothing in quantum PSP has a pre-registered,
sealed, cluster-disjoint benchmark with paired bootstrap CIs and null-calibrated concentration
checks.

### 5.2 The four Sprint 14 results the coordinator asked about specifically

#### A. The best-of-N-from-untrained-circuit control. **SURVIVES. This is the strongest one.**

Searched under several formulations (untrained-circuit baselines, random-initialisation
controls, best-of-N sampling baselines for VQAs, VQA control-plane baselines). **I found no
paper that compares an optimised variational circuit against best-of-the-same-number-of-shots
drawn from the SAME circuit at its initial parameters.** The nearest is Boulebnane et al.'s
uniform-random-sampling control, which is a **weaker** control because it does not hold the
ansatz's induced distribution fixed.

**Why it is the strongest claim we hold.** It is (i) cheap, (ii) universally applicable to every
VQA paper, (iii) a control the field visibly does not run, and (iv) the one that produced our
sharpest result — **0/12, +0.65 to +1.32 A, including from ORACLE warm starts, with the
optimiser erasing its own initialisation.** The combination of "a control nobody runs" and "the
control wins" is a publishable methodological contribution independent of protein folding.
**SURVIVES. Lead with it.** (LITERATURE-SUPPORTED as an absence, with the caveat that absence
searches are weaker than presence searches.)

#### B. The filtering / ordering decomposition of selection. **SURVIVES (RE-SCOPED).**

`sel = pool + FILTERING + ORDERING`, with Legacy at filtering -0.124, ordering **+0.263**, net
+0.139 — and the ORDERING term's sign being "optimise harder, get worse", measured directly.

**Re-tested.** Two adjacent literatures:
* **Docking pipelines routinely separate a "filtering" stage from a "rescoring/ranking" stage**
  as *pipeline architecture* (energy-term filters, then a combined free-energy rescoring).
  That is a workflow distinction, not a decomposition of an outcome.
* **Bias-variance decomposition for ranking** exists in the IR/ML literature, and
  **arXiv:2605.11764 (2026)** performs a variance attribution of a generalisation gap for PROTAC
  activity — and notes that such gaps have been **"documented but not decomposed."**

**Nothing decomposes a SELECTED STRUCTURAL ACCURACY into an additive filtering term and an
ordering term, with signs, against a pool baseline.** **SURVIVES**, but the claim must be stated
as *"an additive decomposition of selected RMSD into filtering and ordering contributions"* and
must not imply we invented the filtering-vs-ranking distinction, which is decades old.
The genuinely novel part is **the ordering term being positive** — i.e. that ordering skill
actively *costs* accuracy — which is the direct measurement of "optimise harder, get worse."

#### C. Tail-restricted discrimination analysis. **PARTIALLY TAKEN. RE-SCOPE HARD.**

> **The concept is TAKEN by: Truchon J-F, Bayly CI, "Evaluating Virtual Screening Methods: Good
> and Bad Metrics for the Early Recognition Problem", *JCIM* 47:488 (2007)** (BEDROC / RIE),
> and the CROC follow-up (*Bioinformatics* 26:1348, 2010).

Virtual screening has had a mature, named theory of **top-of-list vs bulk discrimination** for
nineteen years. BEDROC with `alpha=20` places 80% of the weight on the top 8% of the ranked
list, precisely because ROC/AUC/average-rank are the wrong metrics when only the head matters.
**"Bulk accuracy does not predict top-of-list accuracy" is not a new idea and we must not
present it as one.**

**What survives, and it is specific:**
* The finding that **every objective tested is at chance or worse INSIDE its own low-energy
  tail** (0.512 for the prior; the best *bulk* ranker anti-ranks at **0.390** in its own lowest
  0.1%), and that **the argmin sits at the 41st-48th percentile of its own tail.**
* The methodological point that this must be read against a **random-tail null at 0.524-0.527,
  never 0.500**, and that a gap-matched/binned variant returns to a true 0.500 — a null
  subtlety I found nowhere.
* The three-way agreement `tail mean - tail best = 1.93 A` vs a certified selection gap of
  1.876 and a sampled 2.040.

**Verdict: the FRAMEWORK is taken (cite Truchon & Bayly); the MEASUREMENTS and the tail-null
calibration survive.** Presenting tail-restricted analysis as our invention would be the single
easiest thing for a hostile referee to destroy. **Cite BEDROC in our own Methods.**

#### D. Within-target vs cross-target transfer collapse for a structural objective. **SURVIVES (RE-SCOPED).**

In-band ordering learnable to **0.986 within a target** (overfitting gap -0.0005
[-0.0011,+0.0000]), collapsing to **0.600 across targets** (transfer gap -0.386
[-0.437,-0.332], 770x larger); a **linear** pair potential saturates the within-target problem;
per-target skill correlates **+0.909** with native radius of gyration; **the learning curve
declines.**

**Re-tested.** The random-split-vs-leave-one-target-out collapse is **well documented for
protein-LIGAND scoring and activity prediction** (C17): monotone degradation Random-CV >
Seq-CV > Pfam-CV across 11 MLSFs (*JCIM* 2023); "random splits reward within-target
interpolation, whereas leave-one-target-out measures the novel-target prediction"; and
arXiv:2605.11764 (2026) decomposes such a gap for PROTAC activity, stating that these gaps are
"documented but not decomposed."

**So the phenomenon class is taken. What survives is everything specific:**
* It is measured for a **structural conformational objective** (ordering conformers of one
  peptide), not an affinity/activity regression. I found no instance of that.
* **A linear model saturates the within-target problem**, so capacity, nonlinearity and
  equivariance are all ruled out as remedies — a much stronger statement than "there is a gap."
* **The learning curve DECLINES** — one or two training targets beat twelve. This is the
  opposite of the usual data-hungry story and I found nothing like it.
* **The per-target sign is identified with a physical quantity** (native radius of gyration,
  rho=+0.909), making the failure *explained* rather than merely observed.

**SURVIVES as "a transfer collapse for a structural conformational objective, with the
saturating-linear result, the declining learning curve, and a physical identification of the
per-target sign." Cite C17 for the general phenomenon.**

### 5.3 What is genuinely open, ranked by how defensible it is

1. **The untrained-circuit best-of-N control, and its result.** (§5.2A.) Cheap, universal,
   nobody runs it, and it wins. **Strongest.**
2. **The CVaR estimator/gradient defect triple** — the closed-form tail-baseline bias, the
   non-integer-`alpha*N` upward bias, and the dead-start condition `p(argmin E) >= alpha`.
   Barkoutsos et al. introduce CVaR and never analyse it this way. (§5.1 claim 4.)
3. **The certified global optimum of an ALL-ATOM FORCE FIELD, off lattice, with a CI and a
   filtering/ordering decomposition** — now that Roget et al. own the lattice/contact version.
   (§5.1 claim 3.)
4. **The declining learning curve and the linear-saturation result** for a structural in-band
   objective. (§5.2D.)
5. **The exact torsion-space locality theorem** and its 2-local-impossibility corollary —
   subject to one unrun search pass in the robotics/kinematics literature. (§5.1 claim 1.)
6. **The causal VQE-vs-classical control at matched budget on a molecular objective.**
   (§5.1 claim 5.)
7. **A methodological note naming oracle-selected reporting and proposing
   `pool_best / selected / returned` as a mandatory triple.** (§2.3.) Nearly free.
8. **NEW, and suggested by Roget et al.:** whether the energy-ranking failure is a **length-regime
   artefact of potentials fitted on globular proteins.** Roget's MJ explanation cannot cover
   AMBER ff14SB, which is not knowledge-based. **If AMBER anti-ranks at 9-16 residues while MJ
   anti-ranks for a documented parameterisation reason, the two failures have different causes
   and that contrast is a paper.** HYPOTHESIS.

### 5.4 What must be REMOVED from the story

* **"Nobody has computed the certified global optimum of a folding objective and found it worse
  than random."** FALSE as of June 2026. Roget et al., arXiv:2606.21241.
* **"In-loop CVaR for peptide folding is unclaimed."** FALSE. Uttarkar et al., *PLoS One*
  21(2):e0342012 (2026), plus two 2024 papers from the same group, plus Robert 2021, Kannan 2025.
* **"Tail-restricted discrimination analysis is new."** FALSE. Truchon & Bayly, *JCIM* 47:488
  (2007). Cite BEDROC.
* **"Bulk accuracy not predicting argmin quality is a new observation."** Same citation.
* **"Decoy-bank ranking skill not transferring is our finding."** Handl, Knowles & Lovell,
  *Bioinformatics* 25:1271 (2009). It is a replication.
* **"Spectrum-to-gradient-variance is a new theoretical device."** Fourier/Walsh analysis of VQA
  losses is established (V10) and Kang (2026) gives an exact term-resolved variance
  decomposition. Ours is an application to a molecular energy.
* **Everything S14 already removed** — "first off-lattice quantum folding", "torsion-constrained
  VQE as a novel object", "CVaR as barren-plateau mitigation", "first to notice energy does not
  rank structure". Those removals stand and are reinforced.
* ~~**RESOLVED IN OUR FAVOUR, and the earlier draft of this bullet was wrong.** I initially feared
  that the AMBER/implicit-solvent peptide literature documented ff14SB+GB failing to make the
  native globally stable, which would have made "AMBER anti-ranks the native" a restatement of a
  known defect. **Maffucci & Contini (*JCTC* 12:714, 2016) tests ff96 and the ff99SB series, NOT
  ff14SB** (§2.4). The threat is **downgraded, not eliminated** — cite it as prior context and
  state precisely what ff14SB/GBn2 + enumeration + certified optimum adds. **The residual risk
  sits entirely in the unfetched *PCCP* 2018 study (C4b); fetch it before drafting.**~~

  **⚠ NOT RESOLVED IN OUR FAVOUR. Sprint 16, 2026-09-06 (RETRACT).** The residual risk this bullet
  correctly localised in C4b has **materialised**. Shao & Zhu (*PCCP* 20:7206, 2018) tests
  **FF14SB and FF14SBonlysc with GB-Neck2** — our exact pair — on peptides including a **17-residue**
  β-hairpin, and finds the folding thermodynamics discrepant. The original fear was right.
  **"AMBER anti-ranks the native" is documented expected behaviour of this force field in this
  regime and may not be presented as a discovery.** Re-scoped position and the hostile assessment of
  that re-scoping: §2.4 correction box. The paper is, however, a *mixed* result (it endorses
  FF14SBonlysc/GB-Neck2 as "reasonably balanced" and reports the α-helical 20-mer as working), and
  it reports **no ranking statistic of any kind**, so the observable, the statistics, the AMOEBA
  extension and the stereochemical-repair result all survive.

### 5.5 Claim tier summary

| S14 claim | verdict | authority |
|---|---|---|
| 1. exact locality theorem | **SURVIVES** (one search pass unrun) | negative search |
| 2. spectrum -> gradient variance, no free parameters | **SURVIVES (RE-SCOPED)** — application, not device | V9, V10 |
| 3. certified global optimum on an enumerable space | **TAKEN**; survives only re-scoped to all-atom + off-lattice + CI + decomposition. **⚠ FURTHER NARROWED, Sprint 16 2026-09-06 (RETRACT):** the *force-field* half of the claim is also taken — **Shao & Zhu, *PCCP* 20:7206 (2018)** tests **ff14SB/GB-Neck2** on peptides down to 17 residues and finds the folding thermodynamics discrepant, so "AMBER anti-ranks peptide natives" is documented expected behaviour. Survives only as: the **observable** (rank position in a fixed pool against a certified optimum, which no cited source measures), the **statistics**, the **AMOEBA extension**, and the **stereochemical-repair** result | **Roget et al. arXiv:2606.21241**; **Shao & Zhu *PCCP* 20:7206 (2018)** |
| 4. in-loop CVaR in torsion space | **SURVIVES VERY NARROWLY**; the defect triple is the durable part | **Uttarkar et al. *PLoS One* 21(2):e0342012** |
| 5. causal VQE-vs-classical control at matched budget | **SURVIVES** | negative search; Q7/Q20 are near misses |
| 6. protected statistical instrument | **SURVIVES** | Q1-Q5, X1 all lack it; C18 supports the design |
| A. untrained-circuit best-of-N control | **SURVIVES** — strongest | negative search; Q7 is the near miss |
| B. filtering/ordering decomposition | **SURVIVES (RE-SCOPED)** | docking-pipeline usage is a workflow, not a decomposition |
| C. tail-restricted discrimination | **PARTIALLY TAKEN** — framework taken, measurements survive | **Truchon & Bayly *JCIM* 47:488 (2007)** |
| D. within/cross-target transfer collapse | **SURVIVES (RE-SCOPED)** to a structural objective | C17 owns the general phenomenon |

---

## 6. VERIFICATION STATUS FOR THIS SPRINT

### Primary sources fetched and read in full this sprint
* **arXiv:2606.21241 (Roget, Damour, Cadet, Wang).** PDF fetched, `pdftotext -layout`
  extraction, complete body / Methods / Results / Conclusion read. Every quotation in §1.1,
  §2.1 N1, §2.3 and §3.1 is verbatim from that text. The dataset counts (12,446 / 4,866),
  `d_max=7.8`, `k_min=5`, the enumeration cutoff at length 15, the 10,000-sample Pivot
  procedure, the ~0.05 rho CI amplitude, the ~1.5 A and ~5 A figures, the length-dependence
  result and the HP sensitivity check are all LITERATURE-SUPPORTED.
  **Not verified:** their Figures 4 and 5 numeric values (I read the prose describing them, not
  the figures); the ~1.5 A and ~5 A are their own prose approximations ("around 1.5", "around
  5"), not tabulated numbers.

### Primary sources fetched at full-text or full-abstract level
* **arXiv:2510.06413** (Zhang et al.) — HTML full text. Authors, 75 fragments / 10-14 res /
  PDBbind, 5 candidates each, `E_fuse` formula, the selection sentence, **4.89 mean / 4.70
  median / 1.10 sd**, NetSurfP-3.0, MJ, no CVaR: all confirmed.
* **arXiv:2606.01611** (Yun et al., CD-QAOA) — HTML. APRLRFY, 7 residues, tetrahedral lattice,
  12 qubits + 2 interaction qubits, COBYLA, 10 ns MD, **RMSD 1.39 / 1.90 / 1.70 A vs HF / DFT /
  MD**, dominant bitstring, no CVaR: confirmed.
* **PLoS One 21(2):e0342012** (QuPepFold, Uttarkar et al.) — PMC full text. Tetrahedral lattice,
  `2(N-3)` turn qubits, in-loop `CVaR_alpha`, alpha 0.025/0.05, MJ + geometric penalties,
  1,224 sequences 6-10 aa, ~1.76 M conformers, ~30% faster than expectation-VQE, **no RMSD
  reported**, no classical control: confirmed.
* **bioRxiv 2022.02.17.480937** (Gulsevin & Meiler) — full text. 155 peptides 16-60 aa, the
  per-class RMSD table, both verbatim quotations, the 13/26/18/24/19 rank distribution:
  confirmed.
* **Bioinformatics 25:1271** (Handl et al.) — article page. The 139/149 sentence, the
  three-bias sentence and the rank/z-score sentence returned verbatim: confirmed.
* **Protein Science 10:1470** (Carugo & Pongor) — article page. The length-dependence quotation
  confirmed; **the `rmsd_100` formula is NOT confirmed** (inconsistent forms returned).
* **Quantum 4:256** (Barkoutsos et al.) — journal page. Citation and the CVaR-as-aggregation
  claim confirmed; **alpha values, problem instances and the absence of a random-sampling
  baseline were NOT confirmed from the page.**
* **arXiv:2204.01821** (Boulebnane et al.) — abstract page; abstract quoted verbatim.
  **The quantification of "a small overhead" is UNVERIFIED.**
* **arXiv:2510.15316** (Kannan et al.) — abstract-level. The "minimum RMSD... 1.22 A to 3.11 A"
  wording confirmed this sprint.
* **Maffucci & Contini, *JCTC* 12(2):714-727 (2016)**, doi:10.1021/acs.jctc.5b01211 —
  abstract obtained and quoted verbatim; **full text NOT obtained** (PubMed cookie wall, ACS
  and RSC 403). The critical scope fact — that the tested force fields are **ff96 and the
  ff99SB/ildn/ildn-phi series and NOT ff14SB** — comes from the abstract text itself.

### Search-snippet level only — DO NOT quote as fact without re-fetching
* **C4b, the *PCCP* 2018 AMBER implicit-solvent study (PMID 29480910).** RSC returned 403.
  ScienceDirect returned 403. **Highest-priority unfetched source in this document.**
* The published *Structure* version of Gulsevin & Meiler (588 peptides / 10-40 aa /
  length-normalised CA RMSD). ScienceDirect 403.
* Kang arXiv:2605.01319; the Fourier-expansion VQA literature; Mariyanto *IJQC* 2026; VQEC and
  the augmented-Lagrangian line; Born-machine / CV-DV / Riemannian rows; Truchon & Bayly's
  exact BEDROC algebra; CASP consensus-QA (C16); distance-geometry quotations (C23);
  X2 and X3 (Uttarkar et al. 2024); Q12 and Q21.
* All rows carried from S14 or S13 marked [S] or [P13].

### My own inferences, explicitly not the sources' claims
* **HYPOTHESIS:** that oracle-selected reporting is "endemic and unnamed". Built from five
  fetched papers plus an unsuccessful targeted search. It is an absence claim.
* **HYPOTHESIS:** that AlphaFold2's 13% rank-0 rate is at/below a 20% uniform null. The
  percentages are theirs; the null is my arithmetic.
* **HYPOTHESIS:** that Roget et al.'s MJ-parameterisation explanation cannot cover AMBER
  ff14SB, so a shared anti-ranking would have two different causes. Testable, unrun.
* **HYPOTHESIS:** every "SURVIVES" verdict resting on a negative search. Absence of a search
  hit is weaker evidence than a positive find, and I say so at each claim.

### What I did NOT establish
* I did not fetch C4 (the decisive AMBER/implicit-solvent peptide source).
* I did not fetch Kang arXiv:2605.01319, needed to finalise claim 2.
* I did not read the Barkoutsos CVaR PDF, so the absence of a random baseline there is UNVERIFIED.
* I did not search the robotics / molecular-kinematics literature for the locality result.
* I did not fetch QFold's, Marchand's or Mato's primary texts (still [S]/[P13] from S13).
* I did not obtain X2 or X3 beyond snippets.
* I ran no repository code and touched no benchmark artefact.

---

## STATUS: COMPLETE.
