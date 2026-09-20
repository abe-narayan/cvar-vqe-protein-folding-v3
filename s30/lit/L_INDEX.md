# S30 LANE L -- LITERATURE INDEX (one line per source: verdict, and the information it adds)

Format follows `s29/lit/L_INDEX.md`. The instrument every check is made against: 126 peptides of
9-16 residues, a 500-member retrieval pool whose average carries the statistical content of
**~1.4 independent members** (S30-L7), an ESM-2 650M PCA-32 17-bin distogram trained
leave-fold-out, built chain as the reporting basis, CPU only, 16 GB.

**The standing test applied to every entry:** *what information or mathematical structure does this
contribute that the current pipeline does not already contain?* "Formulated more elegantly" is not
a reason.

S29's carry-forward five (Neyman & Scott 1948; Blau & Michaeli 2018; Brown/Wyatt/Tino 2005 eqs
9-10; McDonald 2023; Cerezo 2025) are **read and not re-read here** -- they are in
`s29/lit/L_INDEX.md` and I build on them rather than repeating them.

---

## Topic 1 -- the reference-state problem at peptide length (`L30_1_reference_state.md`, ledger S30-L6)

| source | information it adds | verdict |
|---|---|---|
| **Shen M, Sali A. Protein Sci 15:2507-2524 (2006)** -- DOPE | the reference state as a *uniform ball of radius a = sqrt(5/3) Rg*, i.e. an explicit, size-parameterised, closed-form `P_ref`; plus the authors' own ablation showing a FIXED reference sphere (DOPE-24) performs substantially worse than the size-adaptive one | **KEPT -- the topic's central source.** Supplies the algebra that makes a statistical potential's compactness loading *forced* rather than measured, and the closed form that prices it |
| Zhou H, Zhou Y. Protein Sci 11:2714 (2002) -- DFIRE, `r^1.61` finite-size reference | the acknowledgement that the ideal-gas exponent is wrong for finite systems | **REJECTED as an import** -- 1.61 was *fitted* to protein-size spheres, so at 9-16 aa it carries the artefact rather than fixing it. NOTED as corroboration |
| GOAP, and "use molecular volume instead of Rg" variants | same mechanism, different size proxy | REJECTED -- no additional information, and volume is worse-conditioned than Rg at this length |
| "a reference state that is too small results in an erroneous preference for loosely packed structures" (survey statement on reference-state construction) | the failure DIRECTION: the reference's size parameter is the knob that sets a potential's compactness preference | KEPT as corroboration -- the knob exists, is published, and is untuned in our channels |
| Yu Z et al. ANDIS. Bioinformatics 35:1499 (2019) | *"native recognition and decoy discrimination cannot be optimized simultaneously with the same parameter sets"* | KEPT (already in the S29 index) -- here it names the PRICE of removing the size term: real between-band signal |
| Fine-grained statistical torsion-angle potentials (137 Ramachandran sub-regions x 20 aa); Ting & Dunbrack (2010) neighbour-dependent Ramachandran | discrimination from **sequence-conditioned** backbone torsion | **REJECTED** -- the conditioner is sequence, measured dead here: memory `phi-carries-no-sequence-signal`, 36.1 deg with full sequence context vs 36.4 deg sequence-blind |
| omega torsion as an additional discriminator | omega's slight (phi,psi) dependence | REJECTED -- omega is fixed at 180 deg by our builder; no information |
| Ca pseudo-torsion / pseudo-angle (theta, tau) potentials, reported to beat DFIRE/dDFIRE/RWPlus on model selection | a scale-free local-geometry channel | **REJECTED AS ALREADY BUILT** -- `s27/ham_lib.py` CAGEO *is* this ("Levitt-style CA-trace potential"), priced at +0.216 partialled in-band by S29-L50. Caught by reading the source before proposing it; this would otherwise have been a re-import, the failure S29-L17 recorded |

Topic 1 count: 8 entries, 3 KEPT, 5 REJECTED. **No importable operator.** The deliverable is a
derivation (S30-L6) and one warning for lane R.

## Topic 2 -- what is identifiable from a pool (`L30_2_common_mode.md`, ledger S30-L7)

| source | what it adds / what it needs that we lack | verdict |
|---|---|---|
| **Kennedy MC, O'Hagan A. JRSS-B 63:425 (2001)** -- calibration with a model-discrepancy term | the canonical formulation in which a systematic shared bias is a first-class object rather than noise | KEPT as framing -- it names our object |
| **Brynjarsdottir J, O'Hagan A. Inverse Problems 30:114007 (2014)** | the calibration parameter and the discrepancy are **not jointly identifiable**; a discrepancy term improves the physical parameter **only** given a strongly informative prior on its SHAPE; and a *wrong* discrepancy prior is worse than none | **KEPT -- the topic's central result.** Proves escape E1's requirement is not optional, and arrives independently at this project's "confidently wrong costs 2-3x absent" |
| bagging / consensus QA lineage (Pcons, ModFOLDclust, DAVIS-EMAconsensus, PWCom) | nothing | REJECTED at family level -- assumes i.i.d. mean-zero member errors with gain `~1/K`; we measure `m_eff = 1.4` out of 75 and `mu != 0` |
| factor models / PCA / ICA on the members | recovers the shared DIRECTION up to the standard sign-and-scale gauge | REJECTED, **and it explains an exact zero**: centring removes `mu`'s component before the decomposition runs, which is why S29-L47's ORACLE global `eta` for the PC1 family is `+0.0000 exactly` |
| blind source separation / ICA identifiability | needs a structural assumption (non-Gaussianity, sparsity, non-negativity) to fix the gauge | REJECTED -- we have none that pins the offset |
| Krogh-Vedelsby ambiguity; Ueda-Nakano bias-variance-covariance (Brown, Wyatt & Tino 2005 eqs 9-10) | **nothing -- it is an identity and it is ours** | KEPT (carried from S29) as the correct framing |
| negative correlation learning; control variates; multifidelity Monte Carlo; Richardson-style two-source contrast | the four inputs S29-L4 enumerated | REJECTED (carried) -- and S30-L7 supplies the *reason* the list is exactly four: each breaks the identification invariance at a different point |

Topic 2 count: 7 entries, 3 KEPT, 4 REJECTED. **The deliverable is a theorem (S30-L7)** plus the
`m_eff = 1.4` conversion factor and a pre-registrable falsifier for lane X.

## Topic 3 -- set functions, subset selection, CVaR on a set (`L30_3_set_selection_cvar.md`, ledger S30-L8)

| source | what it adds | verdict |
|---|---|---|
| Nemhauser, Wolsey & Fisher (1978) `1-1/e`; **Das & Kempe (2011)** submodularity ratio, `1-e^(-gamma)`; randomised greedy for weakly submodular functions | approximation guarantees for cardinality-constrained maximisation | **REJECTED at family level -- every guarantee requires MONOTONE**, and S29-L25 states V is *"neither additive nor monotone"*. Monotonicity, not submodularity, is what we fail first |
| Buchbinder et al. double-greedy (unconstrained non-monotone submodular) | 1/2 without monotonicity | REJECTED -- V is not submodular either; it factors through a centroid |
| determinantal point processes; facility location; diversity-aware selection | a preference over sets (repulsion, coverage) | REJECTED -- V is constant on centroid-equivalence classes, so no diversity structure exists to exploit |
| densest-k-subgraph hardness (Bhaskara et al.); maximum-density-subgraph tractability by max-flow (Goldberg 1984) | the right structural analogy: **fixed-m is the hard formulation, free-m is the easy one** | NOTED -- real, and moot, because the continuous relaxation is cheap |
| Maurey's empirical method / Carathéodory-type approximation | any hull point is within `R/sqrt(m)` of an m-multiset centroid | KEPT -- supplies the relaxation-tightness table (0.22 A at m=75, 1.36 A at m=2) |
| **Maehara T. Oper Res Lett 43:526-529 (2015)** | the CVaR of a stochastic submodular set function **is not submodular**, and admits **no polynomial-time multiplicative approximation** unless P=NP | **KEPT** -- CVaR on a set function changes the complexity class; it is not a decoration |
| **Wilder B. AAAI 2018**; Ohsaka & Yoshida (2017) | the escape is a **portfolio (distribution over sets)**, where `1-1/e` returns via continuous DR-submodular maximisation | KEPT as a design lesson only -- **mismatch stated**: their CVaR is over exogenous randomness, ours over a distribution the optimiser controls |
| Barkoutsos et al. Quantum 4:256 (2020), Prop 5.1 (carried from S29) | CVaR's global-minimiser set `{theta : overlap >= alpha}` is large and flat | KEPT -- and S30-L8 test D adds that **finite shots do not break the degeneracy usefully**: the breaking (+0.016) is smaller than one sd of the shot noise (0.027) |
| the CVaR finite-shot estimator-bias literature (empirical ES / SAA bias) | a plausible mechanism for "objective improves, structure does not move" | **REJECTED BY MY OWN TEST** -- see S30-L8 section 4. Bias is ~0.1% of range at 369 tail shots, does not order states by concentration, does not move the argmin, and runs opposite to my guess |

Topic 3 count: 10 entries, 4 KEPT, 1 NOTED, 5 REJECTED (one of them my own hypothesis). **No
importable operator.** The deliverable is a reframe: *computationally easy, informationally
expensive*.

---

## RUNNING TOTALS

25 sources engaged, 10 KEPT, 1 NOTED, 14 REJECTED across 11 families. **Zero importable
operators**, which is the expected and correct outcome -- the value of this lane is closures and
corrections, not imports.

Three ledger entries: **S30-L6** (reference state), **S30-L7** (common-mode non-identifiability),
**S30-L8** (set selection + a refuted hypothesis of my own).

## THINGS I CHECKED RATHER THAN ASSUMED, BECAUSE THE PROJECT HAS BEEN BURNED FOUR TIMES

- `s27/ham_lib.py` read in source before claiming what our channels' reference states are, and
  before proposing the Ca pseudo-torsion channel -- which turned out to **already exist** as CAGEO.
- `core/quantum.py:179` read before assuming the shot count for the CVaR check.
- The memory note `pool-error-is-68-percent-common-mode` read in its **body**, not its index line.
  The body contains the invariance in prose and the S24 corrections, both of which changed what I
  wrote; the index line alone would have led me to claim the invariance as new.
- `s29/lit/L_INDEX.md` and `s29/s29_L_FINDINGS.md` read in full before any external search, so that
  S29's closures are not re-derived and its near-miss (re-importing PEP-FOLD's retrieval key) is
  not repeated.
- **The two arithmetic scripts this lane cites were moved out of a session scratchpad into
  `s30/lit/` and re-run from their repo paths before the ledger entries were committed.** A cited
  path that exists only in a temp directory is exactly the failure mode `s29_L_FINDINGS.md`
  section 6 and memory `findings-prose-is-not-evidence-of-code` describe.
