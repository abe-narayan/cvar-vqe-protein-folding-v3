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

### 1.5 Part 2, the representation floor on the native's own torsions (ORACLE DIAGNOSTIC; DEMONSTRATED)

`s26/ph_cis.py floor`, `s26/results/ph_cis_floor.json` (complete 126/126), job
`s26/jobs_done/ph_cis_floor.json` (exit 0, 10 s, 0.098 GB). Ledger entry "CIS FLOOR" (see the
ledger tail for its number).

    floor_ca   CA-RMSD, native vs ideal-trans rebuild of its OWN phi/psi      mean 0.347 A (SE 0.029), median 0.272, p90 0.752, max 1.474 (1ID6)
    rebuild_bb the same on N/CA/C (`Peptide.rebuild`)                          mean 0.340 A (SE 0.028), max 1.413
    chain_cost production rmsd_arm - rmsd_avg (built chain minus point cloud) mean 0.166 A (SE 0.018), median 0.098
    Spearman(floor_ca, max omega deviation) +0.828;  Spearman(chain_cost, floor_ca) +0.083;  Spearman(chain_cost, max omega dev) -0.036
    targets with floor_ca > 0.5 A: 32;  floor_ca > chain_cost on 86/126;  cis targets 0 (the cis contrast is empty)

    ST.compare: chain cost minus floor_ca   effect -0.1804  SE 0.0333  MDE 0.0933  iid [-0.2433, -0.1176]  fold [-0.2266, -0.1022]  86W/40L  5/5 folds  (two different objects; "smaller", not "better")

Two readings, both DEMONSTRATED. The constant omega carries a representation cost on this
instrument even with no cis bond: the ideal-trans rebuild of the native's own torsions misses
the native by 0.35 A on average and by more than 0.5 A on a quarter of the targets, and that
miss is omega non-planarity (rho 0.83 with the per-target maximum deviation). And it is NOT
what the projection pays: the 0.166 A chain cost does not correlate with the floor (rho 0.08)
or with the omega deviation (rho -0.04); it is the displacement effect S16 L27 measured, not a
representation effect. Caveat, acted on: `floor_ca` is an UPPER bound on the manifold floor,
because the projection fits phi/psi to a trace rather than rebuilding from the native's
torsions. `s26/PREREG_cis.md` addendum 2 registers `floor2` (the native CA trace projected
through the production projection, `python s26/ph_cis.py floor2`, 10 min CPU) as the tight
number; it runs after the reject jobs and is reported in an addendum here.

---

### 1.6 Part 2b, the tight floor (ORACLE DIAGNOSTIC; DEMONSTRATED; supersedes 1.5's number as "the floor")

`s26/ph_cis.py floor2`, `s26/results/ph_cis_floor2.json` (126/126), ledger entry "CIS FLOOR,
TIGHT FORM". The native CA trace projected through the production projection is 0.083 A from
itself at lam 0.3 (median 0.060, max 0.435, none above 0.5) and 0.043 A at lam 0; the
own-torsion rebuild of 1.5 (0.347 A) overstated the floor by 0.264 A [-0.300, -0.229],
113W/13L, because omega errors accumulate down a rebuilt chain and are absorbed when phi/psi
are fitted to the trace. The constant omega costs about 0.04 A on this instrument; the ramah
prior costs 0.040 A [+0.032, +0.049] on a perfect input (its price where it has nothing to
fix, not a proposal to remove it); the production chain cost exceeds the tight floor by
+0.083 [+0.044, +0.133], so the projection's 0.166 A is the operator's displacement of a
non-native cloud (S16 L27), not representation. Both cis gaps are priced: 0.000 A (no cis
target) and about 0.04 A (non-planarity).

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

### 2.3 How the endpoint will be read, stated before any number exists

The primary reading of "does a steric reject help" comes from arm R (reject and refill to
m = 75): it is the deployable operator, it keeps the averaging mechanism at its production width,
and it cannot empty a set unless the whole 500-candidate pool has no survivor, which happens on
8 targets at 1e4 (20 at 1e3, 1 at 1e5, 0 at 1e6; census, L23). Those targets fall back to the
anchor, exactly as registered, and count as a tie of 0. Arm S (reject, no refill) is the second
reading and it empties 11 sets at 1e4, also to the anchor. Every contrast is then reported
twice: PRIMARY on all 126 (fallbacks included as ties), and as a DECLARED SECONDARY restricted
to the targets whose retained set actually moved (n_reject > 0 and no fallback; at 1e4 that is
116 for R and 113 for S), with n printed beside the MDE. The secondary is a subset chosen by a
native-free property of the operator, not by outcome, and it is not the result; it exists so a
null on all 126 cannot be blamed on the ties. Both were fixed by the coordinator on 2026-09-13
09:10, before `cloud` or `chain` ran.

### 2.4 The point-cloud endpoint (DEMONSTRATED; REFUTED in the harmful direction)

`s26/ph_reject.py cloud` + `report`, `s26/results/ph_reject_cloud.json` (complete 126/126),
`s26/results/ph_reject_report.json`, log `s26/results/ph_reject_report_cloud.log`; ledger L43
(twelve `ST.fmt` blocks verbatim). Basis: POINT CLOUD on both sides; the anchor reproduces
`rmsd_avg` 3.0483 to 1e-6 on every target.

    threshold  arm   vs anchor (all 126)                                    vs matched random (same count)
    1e3        R     +0.532  1.75x MDE  fold [+0.433, +0.657]  36W/69L/21T  WORSE     +0.461  1.87x  WORSE
    1e3        S     +0.196  1.66x MDE  fold [+0.135, +0.233]  29W/71L/26T  WORSE     +0.112  1.20x  WORSE
    1e4  PRIM  R     +0.228  1.17x MDE  fold [+0.147, +0.323]  49W/67L/10T  WORSE*    +0.167  1.26x  WORSE*
    1e4  PRIM  S     +0.108  1.34x MDE  fold [+0.065, +0.153]  38W/75L/13T  WORSE     +0.071  1.00x  NOT MEASURED
    1e5        R     +0.093  0.66x MDE  fold [+0.033, +0.143]  58W/61L/7T   NOT MEAS  +0.042  0.44x  NOT MEASURED
    1e5        S     +0.074  1.17x MDE  fold [+0.039, +0.108]  50W/67L/9T   WORSE*    +0.057  0.95x  NOT MEASURED
    1e6        R     +0.066  0.51x MDE  fold [+0.020, +0.101]  56W/62L/8T   NOT MEAS  +0.030  0.30x  NOT MEASURED
    1e6        S     +0.050  0.97x MDE  fold [+0.026, +0.072]  52W/64L/10T  NOT MEAS  +0.038  0.81x  NOT MEASURED
    (* Type-M zone, 1.0 to 1.3x MDE: sign measured, size an upper bound)
    moved subset at 1e4 (secondary): R n=116 +0.247 WORSE*, R vs RANDR +0.176 WORSE*; S n=113 +0.120 WORSE, S vs RANDS +0.079 at 1.00x
    sweep over the four thresholds (order statistic): R oracle -0.337, split-half -0.147 (44%), k_eff 3.74; S -0.118, -0.056 (47%); the transfer is "the mildest threshold hurts least", not a gain

The falsifier fired the other way. The primary reading (arm R, refill, all 126) is +0.228 A
worse than the shipped top-75 and +0.167 worse than rejecting the same count at random with a
judgment-free refill, 5/5 folds on both; the secondary on the moved targets says the same. The
dose is monotone in the threshold and its limit is doing nothing. The controls locate the harm:
refilling from ranks 76 to 147 costs +0.061 (not measured), shrinking to 35 at random costs
+0.037, and the energy's choice of WHICH members to remove costs the rest. Mechanism from the
census: the condemned members are the ones whose builder-placed side chains clash (96.8%), the
retained set is +0.254 A more expanded in Rg, so the reject keeps expanded members and throws
away compact ones that carried error that cancelled (S23 L5). Power: at 1e5 and 1e6 the arm is
below its MDE (SE 0.05, a gain of 0.13 A would have shown), so those two rungs are UNDERPOWERED
for a small gain and MEASURED against any harm above 0.13 A. No positive result, so no seed
replication; the built-chain run is the cross-basis replication and must agree in sign.

### 2.4b The Adversary's check of L43 (L54, STANDS WITH CAVEAT), answered

All three caveats are accepted and carried forward. (1) The harm is TAIL-CARRIED: at 1e4 the
refill arm's median is +0.003 A against its mean +0.228 (49W/67L/10T, p90 +1.28, worst +4.15 on
8T61); the presentation says "near zero on the median target, catastrophic on the minority
whose pool has no survivor or whose refill reaches deep", never "+0.228 on every target". The
built-chain entry prints the median beside the mean on every line. (2) The +0.228 and +0.167
magnitudes are Type-M-zone (1.17x, 1.26x); "harmful" rests on the measured contrasts (1e3 R
+0.532 at 1.75x, 1e3 S +0.196 at 1.66x, 1e4 S +0.108 at 1.34x, RANDS +0.037 at 1.47x, 5/5 folds
throughout) and on the monotone dose. (3) The "refill is nearly free, the choice costs the
rest" split rests on RANDR vs anchor at 0.57x MDE and is a point-estimate decomposition, not a
measured one; it is downgraded to a reading in section 2.4 and will be re-stated on the built
chain only if RANDR clears its MDE there.

### 2.5 The built-chain endpoint (DEMONSTRATED; the cross-basis replication of 2.4)

`s26/ph_reject.py chain` + `report`, `s26/results/ph_reject_chain.json` (126/126), ledger L86
(eleven `ST.fmt` blocks verbatim), job `ph_reject_chain` (3.3 h, 22.8 projections per target,
0.09 GB). Basis: BUILT CHAIN on both sides; the anchor is the shipped top-75 RE-PROJECTED
through `I.project` (3.2126 A, the leaderboard-rebuild number; the stored production chain is
3.2148 because of its float32 round trip), so every arm shares one instrument; the two chains
differ by more than 0.5 A on 6 targets (2BP4 1.62, 6QAX 1.56, 7JS6 1.54, 2LWU 1.25), the
projection's own branch degeneracy.

    threshold  arm   vs re-projected anchor (all 126, built chain)                          point cloud (2.4)
    1e3        R     +0.563 (median +0.025) 1.77x MDE  fold [+0.441, +0.704]  5/5  WORSE      +0.532
    1e3        S     +0.192 (median +0.001) 1.44x MDE  fold [+0.103, +0.247]  5/5  WORSE      +0.196
    1e4  PRIM  R     +0.248 (median +0.002) 1.17x MDE  fold [+0.120, +0.348]  4/5  WORSE*     +0.228
    1e4  PRIM  S     +0.104 (median +0.003) 1.11x MDE  fold [+0.037, +0.166]  5/5  WORSE*     +0.108
    1e5        R     +0.119 (median  0.000) 0.78x MDE  fold [+0.048, +0.184]  5/5  NOT MEAS   +0.093
    1e6        R     +0.092 (median  0.000) 0.65x MDE  fold [+0.033, +0.138]  4/5  NOT MEAS   +0.066
    R@1e4 vs RANDR +0.157 (1.15x, fold [+0.048, +0.281], 4/5) WORSE*;  S@1e4 vs RANDS +0.055 (0.73x) NOT MEASURED
    moved subset: R +0.270 WORSE*, R vs RANDR +0.169 WORSE*, S +0.116 WORSE*     (* Type-M zone)

Same sign, same monotone dose, same tail structure (medians near zero, p90 +1.4, worst +4.2 on
8T61) as the point cloud; one honest difference, fold 1 flips sign on R@1e4 (-0.002), so R is
4/5 folds here with the fold CI still excluding zero. The projection adds its own noise to
every arm, so the matched-control contrasts that were Type-M on the point cloud are
UNDERPOWERED on the built chain (S vs RANDS 0.73x, R vs PERMR 0.88x); the L54 caveat that the
refill/choice split is a point estimate stands on both bases.

DISPOSITION. The falsifier required R or S at 1e4 to beat the anchor and its matched control on
the built chain; both are worse than the anchor and R is worse than its control. AMBER as a
steric reject at a physical threshold, with or without refill, is CLOSED on both bases in the
two forms the record had not measured, with a harmful sign at 1e3 and 1e4 and an underpowered
null at 1e5 and 1e6. Power on the null rungs: a gain below about 0.14 A on the built chain
would not have been seen; any harm above it would. REFUTED (harmful direction), on both bases.

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

### 3.3 Stage 1, measured (DEMONSTRATED; the toward-member line provisional until replicated)

`s26/ph_c3.py stage1`, `s26/results/ph_c3_stage1.json` (complete 126/126), job
`s26/jobs_done/ph_c3_stage1.json` (exit 0, 10 s, 0.041 GB). Ledger entry "C3 STAGE 1"; the
seven `ST.fmt` blocks are there verbatim and in `s26/C3_RESULT.md`. Bases: arm input the built
chain (`rmsd_arm` 3.2148), arm output the relaxed chain (`rmsd_full` 3.2355); controls displace
the built chain's CA trace by AMBER's own per-atom RMS magnitude (0.220 A).

    AMBER minus do-nothing                      +0.0207  SE 0.0034  MDE 0.0096  fold [+0.0154, +0.0290]  40W/86L  5/5   WORSE   (reproduces production)
    AMBER minus random, same size (16 draws)    +0.0111  SE 0.0035  MDE 0.0099  fold [+0.0062, +0.0171]  52W/74L  5/5   WORSE, Type-M zone (1.12x)
    AMBER minus toward-member, same size        +0.0385  SE 0.0061  MDE 0.0170  fold [+0.0298, +0.0479]  29W/97L  5/5   WORSE
    random minus do-nothing                     +0.0096  SE 0.0015  MDE 0.0042  fold [+0.0077, +0.0122]  29W/97L  5/5   WORSE   (identity predicts +0.0107)
    toward-member minus do-nothing              -0.0178  SE 0.0043  MDE 0.0120  fold [-0.0255, -0.0124]  90W/36L  5/5   BETTER, provisional (replication `ph_c3_stage1_rep`)
    ORACLE cos(AMBER, residual)                 -0.049 (SE 0.015), positive on 36.5%;  random +0.001;  toward-member +0.123
    ORACLE cos: AMBER minus random              -0.0503  MDE 0.0432  fold [-0.0735, -0.0253]  74W/52L   (S16 L27: -0.0491)

The registered prediction (addendum 1: AMBER worse than random by about +0.013, cosine
negative) held. Half of the production step's cost is the size of its move (the random twin
costs +0.0096, the identity's orthogonal prediction +0.0107); the other half is its direction
(cos -0.049, worse than random on 74 of 126 targets). S16 L27 is reproduced on the production
input and operator to the fourth decimal. DECISION RULE: AMBER beats neither matched control,
so "refine with physics" is dropped as an accuracy step and kept as a validity step
(`s26/C3_RESULT.md`).

The one positive line is a control, not the arm: a same-size move toward a random pool member
improves the built chain by 0.018 A (ORACLE cos +0.123). Mechanism on the record: the projection
moved the chain away from the point cloud with a negative cosine (S16 L27) and the members
surround the cloud, so a move back toward any of them recovers part of the 0.166 A projection
cost. It is provisional until the contract's replication (addendum 2: new seeds, reversed
order) lands, and it will be reported as a property of the projection's cost, not of physics.

### 3.4 Stage 2 (ready; waits for lane P's rung file)

`python s26/ph_c3.py probe --input s26/results/p_best_rung_chains.json --pdb <pdb>` under
`--tag AMBER` (est. 1.0 GB), then `stage2` (10 to 30 min, per-target cells). Input rows: `pdb`,
`phi`, `psi` in radians (P's file also carries `ca` and `rung`).

---

---

## 3b. STRAIN DIFFICULTY (tournament rank 4): DEMONSTRATED, REPLICATED; a calibration flag, not a lever

`s26/ph_strain.py`, `s26/results/ph_strain.json` and `_rep.json` (126/126 each), jobs
`ph_strain` / `ph_strain_rep` (35 s, 0.1 GB). Ledger L53 with the full table. Prereg
`s26/PREREG_strain_difficulty.md` (59d8e934, before the run). Basis: `rmsd_arm` (built chain)
as the ORACLE label of the already-emitted structure; every signal native-free and free.

    partial rho(signal, rmsd_arm | n, Rg)   moved +0.433 fold [+0.247, +0.588] 5/5 perm p<0.00025  REPLICATED +0.433 [+0.248, +0.581]
                                            log_e0 +0.241 (4/5 folds)   log_drop +0.250 (4/5)   strain_after +0.133 (CI spans zero)
    rmsd_arm by quartile of moved            2.29 / 2.94 / 3.76 / 3.92 A (32 targets each)
    FAIL18 in the top quartile               null for every signal (best p 0.113)
    not two outliers (without 9KAR, 2BP4: +0.41); not e0 in disguise (moved beats log_e0 by 0.19 in rho)

The first native-free quantity on this record above 0.4 against the per-target error (the
S22/S23 routers and the compactness proxies reached 0.24 to 0.37). It is attachable to every
emitted structure as a confidence flag at zero cost, and it is NOT convertible into a selector
or a weight (the prereg forbids it; S22 L7 / S23 L7 say every such conversion fails held out).
For the presentation: a calibration curve, never a gain.

---

## 3c. BRANCH SELECT (tournament rank 5): REFUTED as an accuracy step; one measured diagnostic

`s26/ph_branch.py`, `s26/results/ph_branch_{solutions,relax,report}.json` (126/126 each; G1 = 0.0
on every target; AMBER jobs `ph_branch_relax` + `ph_branch_relax2`, 0.3 GB, 65 s per target).
Ledger L88 with the five `ST.fmt` blocks. Basis: built chain (`rmsd_arm`) of the chosen
projection solution; anchor = the production choice (re-projected, 3.2126).

    e1 pick minus production          +0.0055  SE 0.011  MDE 0.031  0.18x  39W/45L/42T   NOT MEASURED (a 0.03 A gain would have shown)
    e1 pick minus random branch       -0.1021  SE 0.016  MDE 0.045  2.29x  91W/35L  fold [-0.127, -0.075]  5/5   BETTER
    objective (production) minus random branch   -0.1076  2.65x  5/5           (the objective does the same)
    e0 (raw single point) pick minus production  +0.0922  1.10x  Type-M   WORSE
    ORACLE min over the five: -0.190; valid best-of-5 null 130% of it; split-half transfer 56% (-0.107); k_eff 4.63

The converged relaxed energy chooses among the projection's five branch solutions exactly as
well as the production objective (a 2D torsion prior) and no better: identical pick on 42
targets, a near-equivalent branch on the rest. The unrelaxed single point (the builder's clash,
L23) chooses worse. The branch degeneracy is real and worth about 0.1 A to a perfect chooser;
neither the objective nor the energy is one (each finds the per-target best about 30% of the
time against 20% by chance). CLOSED as an accuracy lever; the report gains the sentence "the
force-field energy, once relaxed, discriminates among the projection's own branches as well as
the torsion prior already does, and the unrelaxed energy does not."

---

## 3d. THE VALIDITY AXIS OF THE PRODUCTION RELAXATION (DEMONSTRATED, native-free)

`s26/ph_validity.py`, `s26/results/ph_validity.json` (126/126; the production relaxation re-run
here reproduces the cache to dCA = 0.0 and dE = 0.0 on every target), ledger L100 (seven
`ST.fmt` blocks), `s26/C3_RESULT.md` addendum 3. Prereg `s26/PREREG_validity_axis.md`.

    axis              built    relaxed   helix     relaxed minus built, fold CI, W/L
    clashes < 2.0 A   0.444    0.008     0.000     -0.437 [-0.555, -0.317]  34W/0L/92T   (34 targets to 1)
    contacts < 2.6 A  3.452    0.119     0.000     -3.333 [-3.992, -2.739]  81W/0L/45T
    min heavy (A)     2.363    2.785     3.079     +0.422 [+0.368, +0.478]  115 of 126 higher
    bond strain       0.000    0.013     0.000     +0.013 (12% on 2BP4)
    angle strain      0.000    0.025     0.000     +0.025
    omega dev (deg)   0.00     6.61      0.00      +6.6 (48 on 1D6X); cis fraction 0.0 to 0.3%
    rama favoured     0.930    0.911     1.000     -0.019, NOT MEASURED (0.84x MDE); outliers +0.008, NOT MEASURED

The relaxation removes the builder's side-chain clashes (the L23 singularity, on the emission
itself) and pays in covalent geometry and peptide-bond planarity, buying no Ramachandran. The
constant helix beats it on every axis, so the claim is the S16 conjunction: clashes removed
while the torsions are held and the trace moves 0.220 A. The presentation sentence is in
`C3_RESULT.md` addendum 3.

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

1. The ensemble prediction was wrong (1.3): 0 of 1,966 models carry a cis bond, not "some".
2. I expected the physical-threshold reject to be surgical. At 1e4 it removes more than half of
   every set, empties 8 whole pools, and is harmful on both bases; the "count set by physics"
   form is "reject most of the set" on this pool.
3. I expected the singularity to mix backbone and side-chain contacts. It is 97% side-chain on
   the condemned members (L23), which made the mandatory reject a test of the builder's rotamer
   placement and moved the question to rotamer_relief.
4. L38's floor (0.347 A) overstated the representation floor by 4x; the registered addendum
   predicted the direction but not the size. The tight floor is 0.083 A (L89).
5. branch_select: I expected null with a harmful sign. The relaxed energy is exactly as good
   as the objective at choosing the branch (beats random by 0.10, 5/5 folds) and no better; the
   raw single point is worse. The mechanism I did not anticipate: relaxation is what makes the
   energy a discriminator among branches, the same fact as S20 L6 from the other side.
6. strain_difficulty: I registered rho 0.2 to 0.3 falling below 0.2 once n and Rg were
   partialled. `moved` came in at +0.43, unchanged by the partial, replicated. The strongest
   native-free correlate of the per-target error on the record, and I had it at 0.3 prior.
7. The C3 toward-member control improves the chain (-0.021, replicated). I expected the two
   controls to be equivalent; a same-size move toward any pool member recovers part of the
   projection's displacement, which is a statement about the projection I had not priced.
8. Two heredocs with lone apostrophes silently failed on the shell (00:46 and 22:45); both
   times the fix was a script file. Recorded because it cost a ledger entry each time.

## 6. WHAT I DID NOT DO AND WHY

- C3 stage 2 (the relaxation on lane P's best C2 rung): waits for `s26/results/
  p_best_rung_chains.json`; the runner (`ph_c3.py probe | stage2`) is ready and tested on the
  production input by stage 1.
- rotamer_relief part B: the probe is queued (`ph_relief_probe`, held at the job cap at 00:30);
  the run follows the probe's peak under the one-AMBER-job rule.
- The reject's built-chain controls carry 4 draws at the primary threshold only (registered);
  16 draws at every threshold would have been 12 h of projection.
- The relaxed emission's heavy-atom panel on lane P's rung output belongs to stage 2.
- No seed replication of the harmful reject results (negatives; the cross-basis run is the
  replication) and none of the descriptive validity panel.

## 7. QUESTIONS FOR THE COORDINATOR

1. L100 (validity axis) and L89 (tight floor) supersede numbers in L24 and L38 respectively;
   both earlier entries stand as written (upper bounds), and I have not retracted them. Say if
   you want a formal supersession entry.
2. The toward-member control's replicated -0.021 A is a property of the projection's
   displacement, not a proposal; if the report wants it as a diagnostic of the projection, the
   wording is in L87.

---

## ADDENDUM A (2026-09-14 02:10) -- RETRACTION IN SECTION 3b PER THE ADVERSARY's L121 (my answer: L123)

The sentence "It is the first native-free quantity in this programme's record with a
correlation above 0.4 to the per-target error of the emitted structure" in section 3b is
RETRACTED, and the framing "the physics reports when the answer is untrustworthy" with it. The
shipped top-75's own pairwise spread, native-free and available before the relaxation, has
partial rho +0.452 with the error (`s26/results/a_strain_vs_spread.json`, fold CI [+0.280,
+0.609], 5/5); `moved` correlates with it at 0.756 and adds +0.082 given it (permutation p
0.39). Restated: the relaxation's displacement is a proxy for the pool's own disagreement,
which is the quantity that predicts the error. The phenomenon, the quartile table and "not a
lever" stand. Section 5 item 6 is amended: the strongest native-free correlate of the error on
the record is the pool spread at +0.45, and the confound list of my prereg should have carried
the pool's own statistics before the operator's (S23 L9, S12).
