# SPRINT 26 -- PHYSICS LANE (PH). FINDINGS.

Status at 2026-09-13 09:00: PHASE 0 NOT SIGNED OFF. Everything below is pre-gate: no RMSD to a
native has been read by this lane. Tiers: DEMONSTRATED / ORACLE DIAGNOSTIC / HYPOTHESIS /
REFUTED / OPEN. Basis is stated at every RMSD. Pre-registrations were written before each run:
`s26/PREREG_amber_reject.md`, `s26/PREREG_cis.md`, `s26/PREREG_c3_control.md` (addenda appended,
never edited).

Artefacts, all `complete: true` on the full key set, provenance-stamped with `ST.save_atomic`:

    s26/results/ph_cis_census.json        126 rows   cis-peptide census (natives, ensembles, pool)   job ph_cis_census     0.038 GB, 90 s
    s26/results/ph_reject_census.json     126 rows   steric-reject census + where the singularity lives   job ph_reject_census  0.117 GB, 100 s
    s26/results/ph_c3_nativefree.json     126 rows   production relaxation displacement + validity scalars   job ph_c3_nativefree  5 s (RSS under-read)
    s26/ph_lib.py  ph_cis.py  ph_reject.py  ph_c3.py  test_ph_synthetic.py   (5 synthetic tests pass; the gate check refuses a native RMSD until the sign-off heading exists)

Ledger: L22 (cis census), L23 (reject census), L24 (C3 native-free).

---

## 0. THE THREE PRE-GATE RESULTS, IN FOUR LINES

1. **There is no cis peptide bond anywhere on the instrument, and it is the database that
   guarantees it.** 0 of 126 model-1 natives, 0 of 1,966 deposited models, 0 of 2,352,893 pool
   windows; the universe's minimum CA-CA step is 3.5045 A, which is `core/data.py:407`'s gate.
   The two-bond-length projection is worth 0.000 A here by construction. DEMONSTRATED (L22).
2. **The steric reject at its primary threshold removes 40 of every 75, empties 11 top-75 sets
   and 8 whole pools, and expands the retained set.** It is not a surgical reject of a few
   impossible members; it is a reject of half the shipped set. DEMONSTRATED, native-free (L23).
3. **96.8% of the candidates AMBER condemns are condemned by a side-chain contact.** Within a
   top-75 the single point is the minimum heavy-atom distance to Spearman -0.74; 40.7 of 75
   rebuilds have an all-atom contact below 2.0 A, against 2.6 of 75 on the backbone plus CB (S19's
   2.66 reproduced). The steric singularity is mostly the deterministic side-chain builder's, not
   the pool's backbones. DEMONSTRATED (L23).
4. **The production relaxation moves the built chain 0.220 A RMS, converges on 125 of 126, and
   breaks a virtual bond on two.** By the S16 identity a move of that size costs 0.0075 A if
   orthogonal to the residual; the step costs +0.0207. Registered prediction for stage 1: AMBER
   worse than a random move of its own size by about +0.013 A. HYPOTHESIS until the gate (L24).

---

## 1. THE CIS-PEPTIDE GAP

### 1.1 The census (ORACLE DIAGNOSTIC on the natives, omega only; native-free on the pool)

`s26/ph_cis.py census`, `s26/results/ph_cis_census.json`. For each of the 126 dev targets, model
1 (the model the instrument scores against), omega_i = dihedral(CA_i, C_i, N_{i+1}, CA_{i+1});
cis iff |omega| < 30 deg; cross-check consecutive CA-CA < 3.3 A; every deposited model
separately; every window of the universe by CA-CA. Sequence check: 126/126 parsed sequences
equal the instrument's.

    model-1 natives with a cis bond (omega)          0 / 126
    model-1 natives with a cis bond (CA-CA)          0 / 126     criteria agree 126/126
    deposited models with a cis bond                 0 / 1,966   (targets with cis in any model: 0)
    universe windows with a step < 3.3 A             0 / 2,352,893
    universe minimum consecutive CA-CA               3.5045 A
    K=500 pools, production top-75 sets with cis     0, 0
    non-planarity |180 - |omega||, 1,507 bonds       mean 1.91 deg, median 0.34, p90 5.77, p99 17.4, max 42.8
    bonds beyond 20 deg / 30 deg                     0.53% / 0.13%   (targets with a bond beyond 30 deg: 1, 9UV5)

### 1.2 Why the answer is zero, with the file and line

`core/data.py:406-407` (`_scan`): a peptide is dropped if any consecutive CA-CA step is below
3.5 A or above 4.1 A. `core/data.py:697-698` (`_extract_fragments`): a fragment window is dropped
if it contains such a step. The 126 dev targets are cluster representatives of `db.load()`
(`s7/debias.py:103-120`), and the sealed benchmark is drawn from the same database, so no target
on either instrument can carry a cis bond, and no retrieved window can either. The measured
universe minimum, 3.5045 A, is the gate's own edge. This is the first time the declared defect
"the projection cannot represent cis" has been priced: on this instrument the price is 0.000 A,
and the reason is the data filter, not the projection.

### 1.3 What damaged my own prediction

I registered that other ensemble models would carry cis bonds on some targets (NMR ensembles are
not filtered model by model). They do not: 0 of 1,966. Either the deposited ensembles of these
short peptides are cis-free or the depositors' models agree on the bond isomer; the census does
not distinguish the two and I do not claim either.

### 1.4 The design note, for the world supply only (no implementation on the production path)

The floor on THIS instrument is the non-planarity tail, priced after the gate (part 2:
`core.data.Peptide.rebuild`, read from `core/data.py:410-419`, is the N/CA/C RMSD of the native
against its own torsions rebuilt at ideal trans geometry; the CA twin is computed beside it).
For a target that does carry a cis bond (the 16 containment-fresh targets, 10 amyloid, memory
`no-fresh-benchmark-exists`), a two-bond-length projection would need:

- Where omega enters. `core/geometry.build_backbone(phi, psi, omega=OMEGA_TRANS)` passes one
  scalar to every `_place_atom` call for CA_{i+1}; `core/project.py` `frames` (line 305) fills
  `t[:, 1::3] = omega` for every bond and `build_ca` / `build_ca_exact` pass the scalar through.
- What a per-bond omega changes. `t[:, 1::3] = omega_vec[None, :]` (length n-1) in `frames`, and
  the same vector in `_place_exact`'s loop in `build_ca_exact`; nothing else. The analytic
  gradient is untouched because omega is not optimised; the bit-exact builder gains one
  broadcast. The CA-CA virtual bond becomes 2.9 A on a cis bond and 3.804 A elsewhere; the
  projection's objective (CA-RMSD to the cloud) is otherwise unchanged.
- Which native-free signal chooses cis. Only the residue after the bond being proline is
  defensible: the literature base rate is about 5% of Xaa-Pro bonds cis and about 0.03% of
  other bonds (Jabs, Weiss and Hilgenfeld, J. Mol. Biol. 286:291, 1999; a literature value, not
  measured here). On this instrument the base rate is 0 of 1,507 bonds and 0 of the proline
  bonds, so no rule can be fitted or validated here; a rule would have to be "cis at Xaa-Pro
  iff the coordinate average's own CA-CA step at that bond is below 3.35 A", which the averaging
  cannot produce from a cis-free pool. The honest statement is that a cis-capable projection
  needs a cis-capable retrieval library first.
- Cost to test. Zero on this instrument (no cis target). On the world supply: a library rebuilt
  with the step gate relaxed for Xaa-Pro bonds, and 16 targets, of which 10 are fibrils, which
  is below any MDE this project can compute. Not proposed.

### 1.5 Part 2 (gated, ready): the representation floor on the native's own torsions

`python s26/ph_cis.py floor` (1 min, CPU). Outputs `s26/results/ph_cis_floor.json`: per target
`rebuild_bb` (N/CA/C), `floor_ca` (CA), `rmsd_arm - rmsd_avg` (the built-chain cost), Spearman
against the maximum omega deviation, and `ST.compare` of the cost against the CA floor. The
cis-vs-non-cis contrast is empty and will not be printed as a subgroup.

---

## 2. AMBER AS A STERIC REJECT FILTER

### 2.1 What the census measured before any RMSD (DEMONSTRATED, native-free)

`s26/ph_reject.py census`, `s26/results/ph_reject_census.json`. On every target the AMBER
cache's pool was asserted equal to `I.pool_idx(u)` bit-for-bit, `amber_verify_max_rel == 0`,
and the production `sub` equal, as a set, to the first 75 of the stable argsort of the cached
distogram score (126/126 pass).

    threshold     pool frac > T   top-75 rejected (mean / median)   zero-reject   all-75 rejected   no survivor in 500   refill depth mean / median   R overlap with anchor   R short of 75
    1e3             0.760           54.5 / 59.5                        1              25                  20                 270 / 234                    0.27                    52
    1e4 PRIMARY     0.586           40.1 / 41.0                        2              11                   8                 201 / 147                    0.46                    14
    1e5             0.446           30.3 / 26.5                        6               3                   1                 167 / 118                    0.60                     6
    1e6             0.345           23.3 / 20.0                        8               2                   0                 135 / 103                    0.69                     2

    retained-set geometry at 1e4, minus the anchor (native-free):   Rg  R +0.254 A (SE 0.058)   S +0.060 (0.012)
                                                                    min |i-j|>=3 CA-CA  R +0.144 (0.028)   S +0.070 (0.013)
    permuted-energy controls reject 44.0 of 75 at 1e4 (vs 40.1 real): the real energy is slightly concentrated in the top-75's tail

The 0.586 is S25's 58.6% (`s25/results/phys_landscape.json`, `AMB_frac_absz_lt_0p1` reads the
same singularity). Targets emptied at 1e4: 1G89 1ID6 1LB7 2MAI 2NB7 2XL1 5MML 5Z5W 7BX2 8UN8
9S5G (top-75); 1G89 1ID6 2MAI 2NB7 2XL1 5MML 7BX2 8UN8 (whole pool). Targets untouched at 1e4:
2MD2, 6A5J. Two of the emptied are in the S16 convergence-gate list (2NB7, 7BX2) and three are
in FAIL18 (1ID6, 2NB7, 1LB7).

Consequences written into the prereg's addendum 1 before any endpoint is read: an empty retained
set falls back to the anchor; the number of fallbacks and the number of targets actually moved
are reported beside every MDE (at 1e4, 124 targets move under R and 113 under S).

### 2.2 Where the singularity lives (DEMONSTRATED, native-free; the optional Part IV measurement)

Every top-75 member was rebuilt with all heavy atoms through the reference builder
(`sidechains.py`, the root module; same construction as `core.amber.build_full_structure`, no
OpenMM), and its closest heavy-atom contact at residue separation >= 2 classified.

    closest contact, all 9,450 members:              bb-bb 1,054   bb-sc 5,300   sc-sc 3,096
    closest contact, the 5,057 members above 1e4:    bb-bb   163   bb-sc 2,627   sc-sc 2,267    (96.8% side-chain-involving)
    members per target with a contact < 2.0 A:       all-atom 40.7 of 75 (SE 1.9)   backbone+CB 2.61 (SE 0.34)   [S19 3.1: 2.66]
    Spearman(e_amber, min heavy-atom distance) in a top-75:   mean -0.743 (SE 0.017), median -0.803
    min heavy-atom distance: above 1e4  1.48 A (SE 0.03);  below 1e4  2.31 A (SE 0.02)

Reading: the AMBER single point on an unrelaxed rebuild is a measurement of the closest atom pair
(rho -0.74), and that pair is a side-chain atom on 97% of the candidates the primary threshold
rejects. The pool's backbones are almost always physically possible (2.6 of 75 backbone+CB
contacts below 2.0 A). So the operator the prereg tests at 1e4 rejects the builder's rotamer
placement. This does not change the falsifier; it fixes, before the endpoint is read, what a
positive result would mean and what a null closes: a null closes "reject the builder's
side-chain clashes", and `s26/IDEA_rotamer_relief.md` part B is the test of the backbone-only
version. Part A of that idea (more than half of the catastrophes side-chain-involving) is
satisfied.

### 2.3 What the endpoint will measure (gated, ready)

`python s26/ph_reject.py cloud` (~5 min): point-cloud RMSD of every arm, every threshold, 16
draws. `python s26/ph_reject.py chain` (~3 h CPU, per-target cells under
`s26/results/ph_reject_chain_cells/`, resumable): built-chain RMSD through `I.project` for the
anchor, R and S at every threshold, and the four controls at 1e4 with 4 draws.
`python s26/ph_reject.py report`: `ST.compare` on all 126 and on the moved subset, fold CI beside
iid, MDE, W/L, concentration null, and the threshold sweep through `ST.best_of_k_within`.
Falsifier and prior: `s26/PREREG_amber_reject.md` sections 4 and 5 (null or harmful).

---

## 3. C3, THE MATCHED-RANDOM CONTROL FOR THE PRODUCTION RELAXATION

### 3.1 The native-free part (DEMONSTRATED)

`s26/ph_c3.py nativefree`, `s26/results/ph_c3_nativefree.json`, from the production cache with
every RMSD key stripped.

    displacement of the built chain by the relaxation, per-atom RMS after superposition   0.220 A (SE 0.008, median 0.197, min 0.103 9WXK, max 0.591 1I93)
    `amber_moved` (restraint RMSD on N/CA/C, unsuperposed)                                  0.233 A
    e0, the built chain's own AMBER energy                                                  median 8.6e4 kcal/mol, min -473, max 1.3e14, 58.7% above 1e4
    e1, after relaxation                                                                    mean -560 (SE 30), max +1262 (9KAR)
    converged, e1 <= 1000 (`core.amber.CONVERGE_MAX_KCAL`)                                  125 / 126
    bond + angle strain after                                                               mean 60 kcal/mol (median 40)
    virtual CA-CA bond, built chain                                                          3.80395 A on every target (sd 1e-16)
    virtual CA-CA bond, relaxed chain                                                        mean 3.867 (SE 0.002); min 3.12 (1M02); max 5.38 (2BP4); 4.86 (9KAR)
    targets with a relaxed CA-CA outside [3.6, 4.0] A                                       6 / 126
    CA contacts < 4.0 A at |i-j| >= 3                                                       built 0.016 per target, relaxed 0.008
    Rg, relaxed minus built                                                                 +0.046 A (SE 0.004)

### 3.2 The registered prediction (HYPOTHESIS until the gate)

S16 section 3.4's identity, n RMSD_after^2 = |r|^2 - 2 v.r + |v|^2, says a move of 0.220 A
orthogonal to the residual costs about 0.220^2 / (2 x 3.215) = 0.0075 A. The production step
costs +0.0207 [+0.0143, +0.0276]. So the matched random control should cost about +0.008 and
AMBER should be worse than it by about +0.013 A, with a negative ORACLE cosine. If stage 1
returns AMBER better than both controls past its MDE with the fold CI excluding zero, the
prediction is wrong and "refine with physics" is an accuracy step; otherwise the presentation
keeps it as a validity step (decision rule, `s26/PREREG_c3_control.md` section 3).

What the validity step can already say with a number: the emitted built chain sits above 1e4
kcal/mol on 58.7% of targets before relaxation (its side chains are placed by the same
deterministic builder that produces the pool's singularity, section 2.2), and the relaxation
brings every target but one below 1000 kcal/mol while moving the CA trace 0.22 A. It also
stretches the virtual bond by 0.06 A on average and breaks it on two targets, which the
heavy-atom validity axis in stage 2 will price properly.

### 3.3 Stage 1 and stage 2 (gated, ready)

`python s26/ph_c3.py stage1` (2 min, CPU): the S16 random displacement (16 draws), the
toward-member displacement (16 draws), AMBER minus each, cosines, on the production chain. Stage
2: `python s26/ph_c3.py probe --input <P's rung json> --pdb <pdb>` under `--tag AMBER`, then
`stage2` (est. 1.0 GB, 10 to 30 min, per-target cells). Input format: rows with `pdb`, `phi`,
`psi` in radians. `s26/C3_RESULT.md` is written when stage 1 is done.

---

## 4. TOURNAMENT IDEAS FILED

`s26/IDEA_amber_reject.md` and `s26/IDEA_cis_peptide.md` (mandatory), `s26/IDEA_branch_select.md`
(the converged relaxed energy picks the projection branch among the multi-start solutions; 1.6 h
AMBER), `s26/IDEA_strain_difficulty.md` (e0, e1, moved as native-free difficulty signals; zero
compute), `s26/IDEA_rotamer_relief.md` (per-residue chi1 relief of the single point; part A
measured and survives; part B 1.2 h AMBER). Each carries a falsifier, an expected effect against
the MDE, a memory estimate and agent-hours.

---

## 5. WHAT DAMAGED MY OWN EXPECTATIONS

1. The ensemble prediction was wrong (section 1.3): 0 of 1,966 models, not "some".
2. I expected the physical-threshold reject to be surgical on many targets (zero rejections on
   clean targets). At 1e4 only two targets are untouched, eight have no survivor in the whole
   pool, and the operator removes more than half of every set. The "count set by physics" form
   is, on this pool, a "reject most of the set" form, and the prereg's addendum says so before
   the endpoint.
3. I expected the singularity to be a mixture of backbone and side-chain contacts. It is 97%
   side-chain on the condemned members. That makes the mandatory reject a test of the builder's
   rotamer placement, and it moves the interesting question to `IDEA_rotamer_relief.md` part B.
4. A shell heredoc with an unbalanced quote silently dropped my ledger entries, STATUS lines and
   prereg addenda at 00:46 and the session ended before I saw it; they were written eight hours
   late (from a script file, not a heredoc). Nothing was lost but time; the census artefacts
   themselves carried their own provenance.

## 6. WHAT I DID NOT DO AND WHY

- No RMSD to a native: the gate is closed. `ph_cis.py floor`, `ph_reject.py cloud` / `chain`,
  `ph_c3.py stage1` refuse to run (`ph_lib.require_gate`) until `PHASE 0 SIGNED OFF` is a ledger
  heading; the check ignores the phrase quoted inside L5's sentence.
- No AMBER compute: the reject uses the 63,000 cached single points; C3 stage 2 waits for lane
  P's rung; `IDEA_branch_select` and `IDEA_rotamer_relief` part B wait for the tournament.
- The minimum heavy-atom distance of the RELAXED emission (the second half of the optional Part
  IV measurement) needs the relaxed backbone, which the production cache does not hold; stage 2
  computes `s16.energy_lib.panel` on input and output and will supply it.
- `IDEA_strain_difficulty` was not run although it costs nothing: it reads `rmsd_arm`, so it
  waits for the gate like everything else.
- The built-chain controls in the reject are 4 draws at the primary threshold only (stated in
  the prereg); if the coordinator prefers 16 draws at every threshold, the chain job becomes ~12 h.

## 7. QUESTIONS FOR THE COORDINATOR

1. The empty-set fallback (anchor) is the rule I registered; if you would rather the emptied
   targets be dropped from the moved-subset contrast only, say so before the endpoint runs.
2. AMBER budget: `IDEA_branch_select` (1.6 h) and `IDEA_rotamer_relief` part B (1.2 h) are only
   worth running if the tournament ranks them; I will not launch either without the ranking.
3. Lane P's rung format for C3 stage 2: I read rows with `pdb`, `phi`, `psi` (radians). If P
   emits CA only, I will project through `I.project` first and say so in the result.
