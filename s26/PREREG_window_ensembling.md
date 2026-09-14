# PREREG_window_ensembling -- TEST-TIME ENSEMBLING OVER RETRIEVAL VARIANTS AT FIXED K (lane P's idea, orphaned to lane W; tournament rank 6; the mandatory test-time-ensembling direction)

Written 2026-09-13 20:05, before any code ran on a real target. Idea file:
`s26/IDEA_window_ensembling.md` (lane P). Code to be written: `s26/w_ensemble.py` (synthetic
tests `s26/w_ensemble_test.py`). Results: `s26/results/w_ensemble_*.json`. Never edited after the
first gated run; addenda appended. Native-free half (`clouds`) stores emissions and reads no
native; gated half (`endpoint`) refuses to run unless "PHASE 0 SIGNED OFF" is in the ledger.

## 0. What the record knows, and exactly how this differs from widening K

- **Widening K (S17 L12, `s17/results/sel_bench.json`)**: ONE shortlist drawn from a wider pool.
  From K = 75 to the full universe the score's top-75 gets better on average (mean 4.287 ->
  3.557) and WORSE at its best (ORACLE best inside the top-75: 2.104 -> 2.572, +0.468 [+0.304,
  +0.642]), because extra plausible windows DISPLACE near-native members out of the one
  shortlist; the consensus readout is best at K = 500 (3.282) and worse at 2000 / full (3.344 /
  3.461). The mechanism is displacement inside a single ranked shortlist.
- **Ensembling at fixed K (this file)**: three shortlists, each of K = 500 and each cut to its
  own top-75 by the shipped score, from three retrieval keys (BLOSUM45 / 62 / 80) or three K
  values (250 / 500 / 1000 at BLOSUM62); each shortlist is averaged in its own medoid frame; the
  three CLOUDS are superposed and averaged; one projection. No shortlist is widened: the
  BLOSUM62 K = 500 top-75 is exactly the production set (its ORACLE best 2.306 A untouched by
  construction), and the other two are independent 75-member shortlists. S17's displacement
  cannot act, because no near-native member is pushed out of any shortlist; what can act is
  averaging three clouds whose errors may or may not be parallel.
- **Why the record predicts null-to-worse anyway, stated before running**: (i) the pool's error
  is 68% common-mode (S23 L9, `s23/results/errdecomp.json` f_common 0.676): a bias shared by
  every member of the universe is shared by all three clouds and cannot cancel; (ii) two
  score-selected candidate sources of DIFFERENT provenance have bias cosine +0.943, above the
  within-source control +0.933 (S24 L3): the shipped score makes the bias parallel whatever the
  source, and three BLOSUM keys over the SAME universe are three sources of the same
  provenance, so the three clouds are predicted parallel and their average is predicted to be
  the shipped cloud plus noise; (iii) the m-ladder is flat on [35, 110] with m = 75 the argmin
  (S12 agg 5a), so a 225-member effective average is predicted slightly worse, not better;
  (iv) retrieval is nearly saturated as a lever: fixing it is worth 0.016 A (S8-6), and the key
  stops mattering by K = 2000 (S7-12). The one route to a gain is BLOSUM45 / 80 retrieving a
  better shortlist on average, bounded by (iv) at about 0.02 A.
- Not previously measured: any combination of two or more score-selected shortlists at fixed K
  (S17 varied K; S12 varied m; S24 L3 compared sources pairwise for bias, never averaged them).

## 1. Hypothesis and exact falsifier

**H_E.** The mean of the three top-75 clouds (BLOSUM45 / 62 / 80 at K = 500), projected once,
is nearer the native on the built chain than the shipped cloud projected once.

**Falsifier.** On the BUILT CHAIN (`arm`, ramah 0.3 multi-start, the 3.2148 basis): `ens3`
minus `shipped` negative beyond its own MDE (`ST.compare`, MDE = 2.8016 x SE), fold-clustered
CI excluding zero, 5/5 folds the same sign, AND `ens3` beats the matched zero-information
ensembling control `boot3` (section 2) by the same standard. Otherwise H_E is refuted at this
instrument with its power stated. The point cloud is carried (stated on both sides); the
selection basis does not exist for an ensemble (no single selected candidate) and is not quoted.
Secondary, reported not decided on: `ensK` (K = 250 / 500 / 1000 at BLOSUM62) under the same
rule; the single-key clouds `b45` and `b80` against `shipped` (are the variant shortlists better
at all); the member-union `union3` (the three top-75s pooled with multiplicity, averaged in the
shipped medoid's frame).

Registered expectation: ens3 minus shipped between 0.00 and +0.03 A on the built chain, i.e.
inside the MDE (about 0.05 A, section 4): the likely outcome is NOT MEASURED with the power
statement "a gain of 0.05 A or more is excluded". Native-free prediction, checkable before the
gate: the three clouds sit within 0.3 A of each other in the median (the shortlists overlap
heavily), so the ensemble cannot move the emission by more than a few tenths on most targets.

## 2. Arms, controls, and the operator, exactly

For each target (universe `s8/generate_univ/<pdb>.npz`: every window's codes `S` in
`core.data.ALPHABET` order, `W` CA coordinates, `sim` = BLOSUM62 sums, `order` = the pinned
stable argsort):

    key k in {45, 62, 80}:  sim_k = M_k[q, S].sum(-1)   (Biopython `substitution_matrices`
                            re-ordered to core.data.ALPHABET; M_62 asserted equal to
                            core.data.BLOSUM62 and sim_62 asserted equal to the universe's `sim`)
    pool_k = core.data.top_k(sim_k, 500)    (the pinned stable argsort; pool_62 == order[:500] asserted)
    top_k  = the 75 lowest shipped Bayes-risk scores in pool_k  (s26/w_selfcopy.emit: production top-75 == cache `sub` asserted for k = 62)
    C_k    = uniform mean of top_k in its own medoid frame       (C_62 == cache `avg_ca` to 1e-13 asserted)

    shipped  = project(C_62)                                             the identity arm (production)
    b45, b80 = project(C_45), project(C_80)                              single variant clouds
    ens3     = project( mean( C_62, superpose(C_45 -> C_62), superpose(C_80 -> C_62) ) )   PRIMARY
    ensK     = the same with C_62 at K = 250 / 500 / 1000 (BLOSUM62)     secondary
    union3   = project( mean over the 225 members (with multiplicity) superposed onto C_62's medoid )   secondary
    boot3    = project( mean of three clouds, each the medoid-frame mean of a 75-member resample
               WITH replacement of the production top-75; seeds s15.seed.stable_rng("w_ensemble", pdb, k) )
               the ZERO-INFORMATION ensembling control matched in cloud count and member count:
               three clouds that carry nothing the shipped cloud does not
    helix    = the constant alpha-helix (phi -63, psi -42 deg), the zero-information scale (S14 L0, 4.0648)

Superposition of clouds: `s12.instrument.superpose_batch` (Kabsch), each variant cloud onto C_62
so the ensemble lives in the production frame. Projection: `s12.instrument.project` (ramah 0.3,
multi-start, `grad="exact"`), one call per arm. Ties: none arise (continuous scores); the top-75
cut uses the stable argsort as production does.

Native-free diagnostics stored before the gate, per target: pool overlap |pool_k & pool_62| / 500,
top-75 overlap, the pairwise RMSD among C_45, C_62, C_80 and among the K-variant clouds, the RMSD
between each arm's chain and the shipped chain (the triangle bound on the change of
RMSD-to-native, L44's tool), and the same for `boot3`.

## 3. Statistics

`ST.compare(arm, shipped, folds=ST.pinned_folds(pdbs), names=pdbs)` for every arm on `arm`
(PRIMARY) and `cloud` (carried); `ST.fmt` verbatim in the ledger; iid and fold CIs; MDE; W/L/T;
median beside mean; concentration null; the verdict from the fold CI and the MDE gate only.
Replication for a positive: a second seed of `boot3` and the reversed fold-processing order
(the ensemble itself has no random element). Power for a null: the MDE of ens3 minus shipped
and the largest gain it excludes.

## 4. Expected effect against the computed MDE

The ensemble differs from the shipped cloud by an averaging of clouds built from the same
universe; the paired sd should resemble a small prior change on the cloud (PREREG_C2 section 4:
MASSFIXW0.1 minus MASS0.0 sd 0.207, MDE 0.052) and the built-chain projection adds its own noise
(rmsd_arm minus rmsd_avg sd 0.205); expected MDE 0.05 to 0.08 A on the built chain. Expected
effect 0.00 to +0.03 A. So the design can refute a gain of about 0.05 A or more and cannot see
a smaller one; the registered deliverable is that power statement plus the native-free cloud
separations.

## 5. Memory, time

Everything comes from the cached universes (one at a time, 1.3 MB each; `s26/w_selfcopy.py
retrieval` ran the same loads and projections at peak 0.116 GB). Per target: three BLOSUM sums
over up to 40,000 windows (milliseconds), three pool scorings, seven projections (shipped, b45,
b80, ens3, ensK, union3, boot3) at ~3.5 s = 25 s; 126 targets ~55 min CPU, checkpointed every
10 targets, est-ram 0.5 GB. One-target probe under jobrun first. Agent-hours: 2 code, 1 run,
1 write-up.

## 6. Operator forks

| fork | DECLARED | NOT TAKEN |
|---|---|---|
| combination | mean of three clouds in the production frame (equal weights) | weights by cloud quality (no native-free quality exists; S22/S23 routers) |
| variant family | three BLOSUM keys at K = 500 (PRIMARY); three K at BLOSUM62 (secondary) | seeds of a random retrieval (there is no random element in retrieval) |
| shortlist rule | the shipped score's top-75 in each pool, m = 75 fixed | a joint top-75 over the union of pools (that is widening K, S17) |
| control | boot3: three resampled clouds of the production top-75 (zero information, matched count) | S25's random-75 null (a different operator: it replaces the shortlist, not the ensemble) |
| basis | built chain PRIMARY, point cloud carried | selection basis (undefined for an ensemble) |
| projection | one projection of the ensemble cloud | projecting each cloud then averaging chains (an average of chains is not on the manifold) |

## ADDENDUM 1 (2026-09-13 22:30) -- what ran; nothing above edited

- Probe `w_ensemble_probe` (1A13; exit 0, 35 s, peak 0.097 GB): gate exact. Run `w_ensemble_clouds`
  (exit 0, 4351 s wall under a four-job load, two governor suspensions, peak 0.113 GB; complete
  126/126, gate 126/126). Endpoint `w_ensemble_endpoint` (exit 0, 10 s, 0.057 GB). Ledger L84;
  findings section 2c.
- Falsifier: H_E refuted; ens3 minus shipped -0.0005 A on the built chain (0.02x MDE 0.0281), ens3
  minus boot3 +0.0042 (0.14x). Power: a gain of 0.028 A or more on the built chain, 0.019 on the
  cloud, is excluded. The native-free prediction (three clouds within 0.3 A of each other in the
  median) held: 0.14 to 0.20 A.
- No deviation. The `sel` basis is not quoted (undefined for an ensemble), as declared.
