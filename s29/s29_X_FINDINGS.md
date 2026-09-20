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

## 4. Results

(filled from `s29/results/s29_X_probe.json` when the probe completes; the ledger entry carries
the ST.fmt blocks verbatim.)

## 5. What damaged my own expectations

(filled with the probe.)

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
