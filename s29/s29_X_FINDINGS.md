# S29 LANE X -- FINDINGS

Standing format: DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS / REFUTED / OPEN, plus "what
damaged my own expectations" and "what I did not do and why". Every number carries its artefact
path. Pre-registration: `s29/PREREG_S29_X.md` (section 0 to 10, addendum 1, addendum 2), all
written before the corresponding numbers existed.

## 1. What this lane built, in one paragraph

A CVaR-VQE whose basis state is not a candidate ID and not a torsion bin, but a CHIMERA: the
chain is cut into S contiguous segments (Rosetta 3-mers, S = min(ceil(n/3), 5)) and each segment
takes its backbone torsions from one of F = 8 retrieved pool members (the DIS top-8 of the
shipped K = 500 pool). Three qubits per segment, q = 3S <= 15, and 2**q = 8**S exactly, so the
register is padding-free and the whole space is exactly enumerable -- every classical control in
this lane is EXACT, not sampled. The Hamiltonian is H = H_diag + H_mix with H_diag the negative
log pseudo-posterior in nats (the shipped distogram's 17-bin posterior consumed as a pair log
score on the built configuration, plus the per-fold Ramachandran log prior as 1-body terms) and
H_mix = -Gamma sum_q X_q a transverse field whose meaning is a transition between chimeras
differing in one bit of one segment's fragment index. The objective is the deployed one plus
that term, F = CVaR_alpha(H_diag; p_theta) - T S(p_theta) - Gamma <psi|sum X_q|psi>, alpha = 0.15,
T = 1 nat, exact statevector, exact parameter-shift gradient on all three terms, Adam, 80
iterations. The readout is the CVaR tail's coordinate average -- R1 uniform (the deployed
operator), R2 tail-mass weighted (the CVaR tail ENSEMBLE, the primary), R3 the state's top-512 --
then the production projection to the built chain.

Code `s29/s29_X_config.py`; tests `tests/test_s29_X.py` (16, green); results
`s29/results/s29_X_probe_<pdb>.json` and `s29/results/s29_X_probe.json`.

## 2. Why this is not a re-run of a closed question (contract rule 10)

The closure table is section 1 of the prereg and is not repeated here. The three lines that
matter: S13 closed the TORSION-BIN space (4 Ramachandran bins per residue) by exhaustive
enumeration under physics energies; S19-S21 closed the BASIN latent (one qubit per residue, von
Mises draws from a prior FITTED to the pool's own marginals) by exhaustive enumeration at n = 126,
where the exact argmin ties a zero-evaluation pool; S28-L21 closed the CANDIDATE-INDEX register
by the set-equality theorem, where the CVaR tail is the classical top-m prefix to 1e-13. The
chimera space is none of those: it is built from the retrieved members themselves, it keeps each
fragment's multi-residue torsion correlations, it CONTAINS its eight parents, and its ORACLE
ceiling had never been measured. Lane T (S29-L15) independently derived that a local mixer -- the
only non-commuting operator class with extensive stable rank -- is meaningful only where basis
states have local structure, which is true here and false in a candidate-index encoding.

## 3. The verification that makes the arithmetic checkable

- A configuration whose S segment choices are all member f rebuilds member f's ideal-geometry CA
  trace with max |difference| exactly 0.0 (`tests/test_s29_X.py::test_parents_reproduce_their_members`).
- H_diag on 3 configurations equals an independent recomputation from the raw `prob` array and the
  raw Ramachandran `cnt` table, written in the test file, to 1e-9
  (`::test_diagonal_hamiltonian_recomputed_independently`).
- At Gamma = 0 the objective and its exact gradient are `core.quantum.free_energy` bit for bit
  (|dF| 0.0e+00, max |dg| 0.0e+00) and the training loop is `core.quantum.run_cvar_vqe` line for
  line (`::test_gamma_zero_is_the_deployed_objective`).
- The transverse field's parameter-shift gradient agrees with central finite differences to
  8.6e-09, and `<psi|sum X_q|psi>` agrees with a direct dense construction to 1e-10.
- The uniform-weight readout is `s12.instrument.coordinate_average` bit for bit (0.0e+00).
- NaN-poisoning `nat_ca` and `oracle_rr` leaves every emitted structure bit-identical.
- Lane T's exact CVaR-optimal law (`cvar_optimal_law`) beats 200 random Dirichlet laws on F and
  reproduces the Boltzmann law exactly at alpha = 1 (TV 0.0), which is T's own assertion.

## 3b. A harness defect in my own code, found by the coordinator (2026-09-20)

Two processes ran the identical unsharded command (`s29X_probe12d` and `s29X_probe12e`;
`s26/jobrun.py` does not deduplicate) and overlapped on one target, 1A13. The artefacts are
undamaged -- all nine parse, all carry exactly 47 arms, 1A13's internals are self-consistent
(q 15, M 32768, segments [3,3,3,3,2], 8 distinct members, 47 unique arm names, both bases on
every arm, the native through the production projection at 0.0282 A) and no stray temp file
exists. But the near-miss is the finding: **`os.replace` is atomic and does not make a SHARED
temp path safe.** Both processes wrote `s29_X_probe_1A13.json.tmp`; an interleave inside that
file would then have been published atomically as a corrupt artefact. Fixed at both write sites
(`tmp = f + ".%d.tmp" % os.getpid()`) and pinned by a test that reads this module's own source,
beside a second test pinning checkpoint-awareness (`os.path.exists(f) and not a.force`), which
is what limited the waste to the single overlapping target. Commit `d6d59487`.

## 4. Results (n = 12, the pre-registered probe; every ST.fmt block verbatim in
## `s29/LEDGER.md :: S29-L56` and `s29/results/s29_X_probe_fmt.txt`)

**DEMONSTRATED (negative), with a mechanism. NO ARM CLEARS ITS BAR.**

1. THE COMPARATOR FIRST (the coordinator's item 1, and it is what makes the ladder readable).
PRODUCTION RESTRICTED TO THESE 12 TARGETS: **3.2529 cloud / 3.3816 chain** against the full-126
anchors 3.0483 / 3.2126, i.e. these targets are +0.2046 / +0.1690 HARDER than the benchmark.
So the ladder is interpretable and the deficit is real rather than a hard-subset artefact.
On the SAME 12 targets (chain / cloud): VQE Gamma_gap R2 3.8388 / 3.7013 - VQE R1 3.8190 /
3.6793 - SA at matched evaluations 3.8574 / 3.8398 - the exact CVaR-optimal law 3.9113 / 3.8912
- the exact argmin of H_diag 3.9122 / 3.9125 - the exact top-75 3.9000 / 3.7939 - the untrained
circuit 3.6937 / 3.6205 - ORACLE best chimera 2.1812 / 2.1801 - ORACLE pool best (K=500)
1.5268. Every native-free arm is 0.31 to 0.57 A ABOVE production on the reporting basis.

2. THE PRE-REGISTERED PRIMARY, P1, AND IT IS THE ONLY ARM ENTITLED TO A VERDICT.

  P1 VQE_g1_s0|R2 - PRODUCTION (chain)
    a 3.8388 (med 3.4080)   b 3.3816 (med 2.7568)   n=12
    effect +0.4572   median +0.4471   SE 0.1388   MDE 0.3890   effect/MDE +1.18
    iid  CI95 [+0.2087, +0.7207]
    fold CI95 [+0.1454, +0.6161]   folds same sign 4/5   per-fold 0:+0.630 1:-0.080 2:+0.165 3:+0.332 4:+0.670
    3W/9L/0T   worst degradation +1.3628 (8HVS)   p90 +1.0789   power 0.91  Type-M 1.05
    VERDICT: WORSE [TYPE-M ZONE: magnitude inflated ~1.05x]

The registered GO rule was "GO to the 126 iff effect <= -0.7 x MDE". It fires in the OPPOSITE
direction: +1.18x MDE, WORSE, fold CI excluding zero, 9 of 12 targets worse. NO 126-TARGET RUN.

3. WHY, AND THE ORACLE SAYS IT WITHOUT ANY SELECTOR: THE SPACE IS POORER THAN ITS PARENT POOL.

  D1p chimera ORACLE best - POOL ORACLE best K=500 (cloud)
    a 2.1801 (med 2.1297)   b 1.5268 (med 1.5795)   n=12
    effect +0.6533   median +0.2468   SE 0.3806   MDE 1.0663   effect/MDE +0.61
    fold CI95 [+0.1130, +1.5173]   folds same sign 5/5   4W/8L/0T   power 0.40  Type-M 1.56
    VERDICT: NOT MEASURED (|effect| 0.6533 <= its own MDE 1.0663, 0.61x)

NOT MEASURED by the MDE rule and 5/5 folds the same sign: recombining the DIS top-8 builds a
space whose BEST structure is 0.65 A worse than the pool's own best on the same targets. The
answer to "is the pool the wrong state space" is that on this construction the pool is the
BETTER one, and no selector can repair a ceiling.

4. D1, WITH LANE D's MATCHED ORDER-STATISTIC NULL (S29-L4 hole (a)), AND THE BRANCH IT FIRES.

  D1  chimera ORACLE best - SCRAMBLED ORACLE best (cloud)   effect -0.1086  MDE 0.1796  -0.60x
      fold CI95 [-0.1559, -0.0419]  5/5 folds  9W/3L   VERDICT: NOT MEASURED
  D1u chimera ORACLE best - best PARENT (cloud, UNMATCHED)  effect -0.8065  MDE 0.7547  -1.07x
      fold CI95 [-1.2435, -0.3830]  5/5 folds  12W/0L   VERDICT: BETTER [TYPE-M 1.09]

**87% of the apparent value of recombination is the order statistic.** Against the eight
parents the space looks worth -0.807 A; against a SCRAMBLED space of identical cardinality,
identical parents and identical marginal fragment content it is worth -0.109 A and NOT MEASURED.
Lane D required that control before any number was read, and it is decisive. Prereg addendum 1
branch (b) therefore fires: recombination is worth less than 0.1 A against the matched null, so
the space is no richer than its parents and nothing further is built on it.

5. LANE T's LADDER (S29-L15), IN ORDER, WITH THE GATE-1 PASS EXPLAINED RATHER THAN CELEBRATED.

GATE 1 (TV(p_Gamma, p_0) > 0.45 on the sampled distribution) CLEARS: mean TV 0.5765 / 0.5631 /
0.5563 at 0.5x / 1x / 2x Gamma_gap, above 0.45 on 11/12, 11/12, 10/12. **But the diagnostic that
decides how to read it says the pass is FLATTENING, not structure**: at Gamma_gap the trained
state's entropy is 95.7% of maximal and <sum_q X_q> is 88.2% of its maximum q, i.e. the state is
close to |+>^q -- which is an EIGENVECTOR of H_mix, exactly the disqualifier lane L's S29-L13
names ("the prepared object must not be an eigenvector"). The mixer at the native-free coupling
scale overwhelms H_diag rather than competing with it. Across the grid the realised tail size
goes m/(alpha M) = 0.0225 -> 0.3792 -> 0.6476 -> 0.7868: **the mixer's whole endpoint channel is
m**, which is lane T's Q1 reduction surviving intact into configuration space. (The tail is an
energy PREFIX and EQUALS the classical top-m on 12/12 at every Gamma including 2x. That is NOT a
discovery: `cvar_from_probs` allocates the alpha mass along the energy order, so a full-support
state's tail is a prefix BY CONSTRUCTION -- lane L, S29-L13, Barkoutsos eq 12. It is reported
because it shows the mixer produced no exact zeros, so m is the only channel it left open.)

GATE 2 (the correctly-named classical counterpart: SA over the same space at 2**q evaluations,
tail-averaged through the identical extraction -- NOT an eigensolver):

  P2 VQE R2 - SA R2 (chain)
    a 3.8388 (med 3.4080)   b 3.8574 (med 3.5385)   n=12
    effect -0.0186   median -0.0350   SE 0.0713   MDE 0.1999   effect/MDE -0.09
    fold CI95 [-0.1023, +0.1179]   folds same sign 2/5   7W/5L/0T   power 0.06  Type-M 9.02
    VERDICT: NOT MEASURED (|effect| 0.0186 <= its own MDE 0.1999, 0.09x)

GATE 3 is section 2 above and it fires the wrong way. The ladder therefore ends here.

6. THE TEN CONTROLS (contract rule 15) AND THE TWO LANE-D ADDITIONS. Built chain, n = 12:

    P3   VQE R2 - VQE R1 (the quantum stage beyond m)          +0.0198  MDE 0.1051  +0.19x
    P3c  VQE R2 - PR-matched random weights (lane D hole (c))  +0.0271  MDE 0.1029  +0.26x
    P4   VQE R2 - GIBBS at matched ENTROPY                     -0.0254  MDE 0.1182  -0.22x
    P5   VQE Gamma - VQE Gamma = 0                             -0.1171  MDE 0.2726  -0.43x
    M6   VQE(Gamma=0) - the EXACT CVaR-optimal law             +0.0446  MDE 0.0894  +0.50x
    C-untrained   VQE - the untrained circuit                  +0.1450  MDE 0.3095  +0.47x
    C-product     VQE - the product-state restriction          -0.0355  MDE 0.1463  -0.24x
    C-seed        seed 0 - seed 1                              +0.0426  MDE 0.1056  +0.40x
    C-diagonalised VQE - the exact ground state                -0.0380  MDE 0.1151  -0.33x
    C-perm        VQE - VQE on a permuted posterior            +0.1319  MDE 0.2761  +0.48x
    C-ordstat     VQE - the exact classical top-m              +0.0228  MDE 0.1057  +0.22x
    BESTOFN       best-of-N(untrained, by F) - VQE             -0.0088  MDE 0.1308  -0.07x
    BESTOFN       best-of-N(untrained, by F) - 1 untrained draw +0.1362 MDE 0.2869  +0.47x

**ALL THIRTEEN ARE NOT MEASURED**, every |effect| between 0.009 and 0.145 A against MDEs of
0.089 to 0.369. The strongest single statement in the lane is BESTOFN, the control project
memory demands (`concentration-is-wrong-when-discrimination-binds`): at matched budget,
**best-of-N from the UNTRAINED circuit selected by F is indistinguishable from the trained VQE
(-0.0088, 0.07x MDE)**, while the optimiser demonstrably works -- F falls by half or more on
every arm (e.g. 117.3 -> 73.9, 88.1 -> 82.0, 325.1 -> 108.5). This RE-TESTS S20's scope fix
rather than assuming it: on this continuous configuration encoding the circuit trains, and the
emitted structure is still not separable from choosing the best of 80 random parameter draws by
the same objective. Nothing the quantum stage does survives to the chain.

7. MECHANISM BESIDE OUTCOME (contract rule 18; addendum 20(c)). H_diag does not order this space
for nativeness: corr(H_diag, RMSD) = +0.0199 over all configurations, the ORACLE-best chimera
sits at the 62nd percentile of H_diag's own ordering, and the native's pair term is WORSE than
71.7% of chimeras (79.3% after the native is put through the production projection, so it is not
a manifold artefact). The exact argmin of H_diag emits 3.9125 where the space's mean is 3.8771 --
**the objective's exact optimum is no better than a random point in the space it just enumerated
in full**, which is S13 §4 and S21 L14 reproduced in a third space. Shape: emitted Rg 6.72 and
mean virtual bond 3.288 for the VQE against 3.804 for the ORACLE best chimera -- the emitted
structures are contracted, as every averaging readout here is. Reachability: the trained state's
mean overlap with the exact ground state of diag(H_diag) - Gamma sum X is 0.4355 (range 0.0003
to 0.8895), so the circuit is far from that state on some targets -- and it does not matter,
because the diagonalised equivalent (C-diagonalised, -0.0380) is itself null.

8. MULTIPLICITY (contract rule 17). This family is 125 ST.fmt blocks: 61 built-chain contrasts,
48 of them arm-vs-production, plus 63 cloud blocks, over 47 arms carrying a chain number. With
K = 48 one-sided contrasts, P(at least one spurious at alpha = 0.05) = 0.9147 and the Bonferroni
two-sided bar is z = 3.279 against 1.96. Only P1 was pre-specified as the primary and it is the
only arm read as a result. The per-target MAXIMUM over the 21 native-free deployable arms is
priced with `ST.best_of_k_within`: observed gain -0.4476 A, **best-of-k null -0.6088 A, share
accounted 136%, k_eff 6.93, split-half transfer -0.0312 A (7.0% of the oracle gain)**. The
per-target minimum is entirely an order statistic and nothing transfers -- exactly the audit
`grid-oracles-are-order-statistics` requires.

9. ARTEFACT INTEGRITY, AND A HARNESS DEFECT IN MY OWN CODE (coordinator's item 4). All 12 files
parse with exactly 47 unique arm names and both bases; `s29_X_probe_1A13.json` (the target two
processes were briefly writing concurrently) is intact -- 47 arms, q 15, M 32768, sha256
07eb9eebfbb95929..., no stray temp files. Nothing is dropped. But the near-miss was real:
**`os.replace` is atomic and does NOT make a SHARED temp path safe.** Both processes wrote
`s29_X_probe_1A13.json.tmp`, and an interleave inside that file would have been published
atomically as a corrupt artefact. Fixed at both write sites (`tmp = f + ".%d.tmp" % os.getpid()`)
with tests pinning it and pinning checkpoint-awareness (which is what limited the duplicate run's
waste to the single overlapping target). Commit d6d59487.

10. THE R3 READING IS REFUTED (coordinator's item 2). On the 5-target read UNTRAINED|R3 was the
best non-oracle arm by 0.17 to 0.22 A. At n = 12 it is +0.1358 / +0.0247 / +0.1749 / +0.1721
(cloud) against the four trained R3 arms, i.e. 0.08x to 0.55x MDE, all NOT MEASURED, with the
seed-1 cell at 6W/6L. And R3 is the LEAST informative readout, not the most: arm-spread divided
by the mean paired sd is R1 1.047, R2 0.753, **R3 0.459**, and R3 captures a mean 0.3116 of the
state's probability mass -- it is a top-512 readout, not a full-state one (lane D's S29-L4 hole
(f), now quantified). It separated at n = 5 because it is the noisiest, which is the failure
mode the coordinator named.

11. THE NINE QUESTIONS (charter section 11), answered with measured values. (1) Basis state: a
chimera, i.e. which of 8 retrieved members supplies each of S segments; 8**S real-fragment
recombinations, the 8 parents among them (verified bit-exact). (2) Hamiltonian: the negative log
pseudo-posterior of the chimera under the shipped distogram plus the Ramachandran log prior;
H_mix is the one-bit transition operator. (3) Why it should correlate with structure: the
distogram orders the bulk (S13, S21 L14) -- measured here, it does not: corr +0.0199, ORACLE-best
at the 62nd percentile. (4) CVaR optimises the 0.15 tail of H_diag under p_theta minus T x
entropy minus Gamma x coherence. (5) The ansatz has 27 to 45 real parameters over 512 to 32,768
amplitudes; its mean overlap with the exact ground state is 0.4355. (6) The optimiser reaches
low F (F halves or better on every arm) and the exact CVaR-optimal law is reached to within
+0.0446 A of emitted RMSD. (7) What the quantum output could contain: a correlated superposition
over compatible chimeras -- measured at 88.2% of the mixer's maximum, i.e. nearly the uniform
product state |+>^q, so it contains no correlation. (8) A classical control reproduces it:
GIBBS at matched entropy, SA at matched evaluations, the exact top-m, the product state, and
best-of-N all tie it (13 nulls). (9) The built chain moves the WRONG way: P1 = +0.4572,
1.18x MDE, WORSE.

VERDICT. **The divergent direction is CLOSED, negatively, with a mechanism.** (a) The chimera
space is poorer than its parent pool at the ORACLE (+0.6533 A, 5/5 folds) and recombination is
worth -0.109 A against its matched order-statistic null, so 87% of its apparent value was
best-of-8**S. (b) The endpoint is +0.4572 A WORSE than production at 1.18x MDE, and the GO rule
fires in the opposite direction: no 126-target run. (c) The local transverse field -- lane T's
smallest qualifying non-commuting operator -- passes the TV gate only by driving the state to
|+>^q, an eigenvector of the mixer, which is lane L's own disqualifier; its entire endpoint
channel is the tail size m, so lane T's Q1 reduction holds in configuration space as well as in
the candidate-index encoding. (d) With the control project memory demands, best-of-N from the
untrained circuit at matched budget is indistinguishable from the trained VQE (0.07x MDE) even
though the optimiser plainly works. NOTHING IN THIS FAMILY CLEARS ITS BAR.

WHAT DAMAGED MY OWN EXPECTATIONS. I registered that recombination would be worth 0.3 to 0.8 A at
the ORACLE over the top-8's best member. Against the parents it measures -0.807 A, which looks
like a hit -- and lane D's matched null, demanded before I read a number, shows -0.109 A. My
prior was right about the arithmetic and wrong about the science, and it would have been a
positive in the record if the control had not been required first. Second: I expected a
non-commuting term to break the tail-is-a-prefix reduction. It does not, for a reason I should
have derived rather than measured -- the extraction is a function of p and the energy ORDER
alone, so non-commutativity enters only through p, and p only moves m.

SCOPE, STATED SO IT IS NOT OVER-READ. This closes F = 8 DIS-top parents x contiguous 3-mer
segments x q <= 15, with H_diag = the pair log score plus the Ramachandran prior, on 12 targets.
It does not close fragment recombination in general (a different parent set, non-contiguous or
overlapping segments, or more parents than 8 are untested), and the 12-target probe is never
quoted as evidence for the instrument (contract rule 16). What it does close is this lane's own
hypothesis, and the reason is a CEILING, not a selector: no operator can recover 0.65 A that the
space does not contain.

COST-RMSD METER (contract rule 19), run by lane D before the endpoint (S29-L6): the pair log
score is the first cost in the record that does NOT anti-order the near-native ladder (+0.200
rho above the shipped cost, fold CI [+0.100, +0.343], 5/5 folds) and it still does not recognise
nativeness (native at the 37.8th percentile of its own pool). Its gradient is undefined on
114/126 targets because a binned log score is a step function -- which does not bind here,
because the search is discrete and exhaustive, and is stated so the number is not reused out of
scope. Contract addendum 21: H_diag is inside the bounded marginal class, so no gradient or
cosine claim is made; what was tested is ORDERING over a different set and NON-CONTRACTION, and
both fail.

ARTEFACTS. `s29/results/s29_X_probe.json` (the aggregate, every ST.fmt block under `fmt`);
`s29/results/s29_X_probe_<pdb>.json` x 12; `s29/results/s29_X_bestofn.json` and
`s29_X_bestofn_<pdb>.json` x 12; `s29/results/s29_X_probe_fmt.txt` (the ladder in lane T's
order); `s29/results/s29_X_bestofn_fmt.txt`; code `s29/s29_X_config.py`, `s29/s29_X_report.py`;
tests `tests/test_s29_X.py` (19, green); pre-registration `s29/PREREG_S29_X.md` with addenda 1
to 3; findings `s29/s29_X_FINDINGS.md`. Jobs `s29X_probe1`, `s29X_probe12e`, `s29X_bestofn`
under `s26/jobs_done/`.

## 5. What damaged my own expectations

1. **I registered that recombination would be worth 0.3 to 0.8 A at the ORACLE over the
   top-8's best member, and against the parents it measures -0.807 A** -- my prior, met.
   Against lane D's matched SCRAMBLED null (identical cardinality, identical parents,
   identical marginal fragment content) it is -0.109 A and NOT MEASURED: 87% of it was
   best-of-8^S. Had that control not been demanded before I read a number (S29-L4 hole (a)),
   this lane would have put a 0.8 A "recombination is worth something" claim into the record.
2. **I expected a non-commuting term to break the tail-is-a-prefix reduction, and it does not**
   -- 12/12 targets at every Gamma including 2x Gamma_gap. I should have DERIVED that rather
   than measured it: the extraction is a function of p and the ENERGY ORDER alone, so
   non-commutativity enters only through p, and p only moves m.
3. **I expected gate 1's TV pass to be the interesting result.** It is the opposite: the pass
   is produced by the mixer flattening the state to ~|+>^q, an EIGENVECTOR of H_mix, which is
   precisely the disqualifier lane L had already derived (S29-L13). My own gate was passed for
   a reason that fails the condition the gate exists to protect.
4. **The 5-target read looked like a finding and was not.** UNTRAINED|R3 as the best non-oracle
   arm collapses to 0.08x-0.55x MDE at n = 12, and R3 is the LEAST informative of the three
   readouts (separation-to-noise 0.459 against R1's 1.047), not the most.

## 6. What I did not do, and why

- The 126-target instrument: contract rule 16 forbids it before the probe, and the probe's GO
  rule is pre-registered.
- AMBER, anywhere: no physics energy ranks nativeness on this pool (S25 L16), and it is not part
  of this hypothesis.
- Encoding (i), the per-residue torsion bins: closed twice by exhaustive enumeration (S13, S21),
  and the brief names (ii) as the one the record has never run.
- Any tuning of alpha, T, Gamma, F, the segment length, the member rule or the readout on an
  RMSD. The one compute-driven choice (the register cap S <= 5) is in prereg addendum 1 with the
  measured cost that forced it.
- A shot-noise study: the statevector is exact at these register sizes, and an exact simulator is
  never a cause (contract rule 9).
