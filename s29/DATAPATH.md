# S29 DATAPATH -- sequence to built chain, stage by stage (lane M, 2026-09-19)

The full path from a target sequence to the built chain that `s27/results/chain_rows.jsonl :: DIS`
scores at 3.2126 A, written from the code, with the file and line of every function, the shape and
unit of every object, what each stage CREATES, TRANSFORMS, COMPRESSES or DESTROYS, and the
dimension or entropy count before and after where it can be stated. Line numbers are of the working
tree at HEAD `b1618a6a` (production code frozen at `a15406c`). Every count marked "measured" was
computed by lane M on 2026-09-19 with the commands in `s29/s29_M_harness_audit.md` or the one-off
script recorded in `s29/s29_M_FINDINGS.md`.

**Read this first.** Two things about the path that the architecture diagram (`ARCHITECTURE.md`
section 1) does not make plain:

1. **The anchor 3.2126 A never passes through the quantum stage.** `s27/run_vqe_chain.py --chain`
   (lines 124-130) takes the top-75 of the tie-safe energy order, averages, projects and scores;
   `arm_vqe` is called only in `--vqe` mode (lines 111-123, `vqe_rows.jsonl`). The production
   `core.pipeline.Config` has `quantum: bool = False` (`core/pipeline.py:179`), and the production
   cache record that produced 3.2148 (`bench_results/cache/1fc9f2dcf489e2fb/1A13.json`) carries
   `quantum: null, q_ca: null, n_top: 75`. The CVaR-VQE is a research arm evaluated beside the
   spine (S25 suite, S27 `vqe_rows.jsonl`, S28 lanes); by the set-equality theorem its emitted set
   is the classical top-m of the same energy (`s24/d_harness.py:32-62`), so the two arms differ
   only by the rung m (74/71 at seeds 0/1 against the fixed 75; S28-L21).
2. **The only place three-dimensional coordinates are CREATED is retrieval (stage 5).** Every
   later stage is a selection, an average, or a projection of real deposited backbone geometry.
   The distogram (stage 3) creates a distribution over distances, never a coordinate.

Notation: n = target length (9 to 16); npairs = (n-1)(n-2)/2 = 28..105 (pairs with j-i >= 2,
`s12/instrument.py:128 pair_index(min_sep=2)`); K = 500; M = 75; nw = number of library windows.

---

## Stage 0 -- the corpus, the identity clusters, the pinned folds, the 126 targets

| item | where |
|---|---|
| peptide bank | `core/data.py:421 _scan`, `457 build_peptides` -> `peptide_db.npz` |
| fragment bank | `core/data.py:715 _extract_fragments`, `746 build_fragments` -> `fragment_db.npz` |
| identity | `core/data.py:182 identity` (Needleman-Wunsch matches / LONGER length; declared leak) |
| clusters | `core/data.py:518 clusters` (single linkage at `IDENTITY_THRESHOLD = 0.6`, line 404) -> `peptide_clusters.json` |
| folds | `core/data.py:558 folds` (cluster ids shuffled with seed 0, i mod 5) -> `peptide_folds.json` (PINNED) |
| fold-safe fragments | `core/data.py:794 fold_fragments` (drop fragments at identity >= 0.6 to any held-out peptide of that fold) |
| the 126 targets | `s7/debias.py:103 tuning_targets` (one peptide per cluster, 9 <= n <= 16, cluster-disjoint from `dev_set(24)` and `benchmark()`), materialised once as `s8/generate_univ/<pdb>.npz` and read back by `s12/instrument.py:41 targets()` in pdb-sorted order |

Inputs: `pdbs/` + `pdbs_ext/` deposits (model 1, `core/geometry.py:823 native_coords_from_pdb`),
`prots/` (1,001 protein deposits). Filters: length 8..26 (`core/data.py:398`), all torsions finite,
CA-CA step in [3.5, 4.1] A, and rebuild-from-own-torsions backbone RMSD <= 1.5 A (`REBUILD_TOL`,
line 403); fragments: lengths 9..20 at stride 5, at most 6 per protein, rebuild <= 1.0 A (lines
707-709, 746).

Objects (measured 2026-09-19): **787 peptides** (8..25 aa) and **6,003 fragments** (9..20 aa); each
entry is `seq`, `ca (n,3)` float64 A, `phi/psi (n,)` radians. Fold sizes 158/140/179/161/149
peptides. `fold_fragments(f)` keeps 5,988..6,003 of 6,003: **the fragment bank is shared by all five
folds** (0..15 fragments removed per fold), so "leave-fold-out" is a statement about the 787
peptides only. The 126 targets: 126 distinct pdb ids, 126 distinct sequences, lengths 9..16, fold
counts 25/23/25/23/30; every target's universe `fold` equals `peptide_folds.json[seq]` (0 mismatches).

CREATED: the sequence-structure pairs the whole system learns from and retrieves from, out of
deposited coordinates. DESTROYED: everything but CA positions and (phi, psi); side chains, the
carbonyl geometry, partner chains, ligands, the ensemble of NMR models beyond model 1 (`model_index=0`,
`core/geometry.py:808`), cis peptides and distorted valence geometry (the rebuild filters). The
identity normalisation by the longer sequence lets a verbatim k-mer inside a longer entry fall below
0.6 and land in another fold (ARCHITECTURE section 2.1: 4/126 dev targets carry a verbatim self-copy in
their own fold model's training set). `tuning_targets` calls `db.benchmark()`, which reads
`results/benchmark_manifest.json` (`core/data.py:588, 655`) -- at DEFINITION time only; the run-time
path reads the pinned universes and never opens the manifest (audit check 3).

Entropy: a target sequence carries n log2 20 = 38.9..69.2 bits; the peptide bank 787 sequences,
the fragment bank 6,003.

---

## Stage 1 -- sequence -> ESM-2 embedding and contact map

`core/data.py:938 esm_compute` runs `esm2_t33_650M_UR50D` (`ESM_MODEL`, line 886), final layer:
per-residue representation **(n, 1280) float32** and the contact head **(n, n)** probabilities
(`esm_contacts`, line 1025 -- "the only part of ESM-2 ever supervised on structure"). Served from
`esm_small.npz` via `esm_raw` (970).

`core/data.py:989 fit_esm_pca`: one GLOBAL PCA over the corpus residues (mean `mu`, loadings
`W (1280, 32)`, `scale`), `N_PCA = 32` (line 887). `1007 esm_embed`: `((rep - mu) @ W) / scale`
-> **(n, 32) whitened**.

COMPRESSED: 1280 -> 32 per residue, a factor of 40, by a fixed projection fitted on the corpus (not
per fold; S26 L67 measured a per-fold refit at +0.018 A built chain, 0.13x MDE, so the global fit is
not a leak worth anything; S26 L72 measured 128 components at +0.076 A, 0.43x MDE). The contact map
is consumed only through five gathers per pair and two row statistics (stage 2).

---

## Stage 2 -- pair features

`core/predict.py:94 pair_features` -> base block, 42 features per pair (measured): 7 separation and
position terms (`s, log s, s/n, n, i/n, j/n, min(i, n-1-j)`, lines 112-114) + 5 physicochemical
properties (`_PROPS`, lines 64-75: Chou-Fasman helix, sheet, Kyte-Doolittle hydropathy, charge,
volume) in 7 combinations (`P[i], P[j], W[i], W[j], P[i]*P[j], |P[i]-P[j]|, between`, line 115).

`core/predict.py:129 features` appends the ESM block, 141 features: `con[i,j]`, its 4 diagonal
neighbours, `deg[i], deg[j]` (row means of the contact map), 3 separation-shell contact masses for
each of i and j, and `emb[i], emb[j], emb[i]*emb[j], |emb[i]-emb[j]|` (4 x 32 = 128) (lines 141-160).
**Total width 183** (`feature_width`, line 169; measured 183). Output: `X (npairs, 183) float32`,
`i, j (npairs,)`.

TRANSFORMED: sequence + embedding -> a per-pair vector. DESTROYED: pairs with j - i = 1 (never
predicted, never scored; `min_sep = 2`, `core/pipeline.py:133`); every residue's context beyond +-2
(the `W` window, lines 100-101) except through the ESM vectors.

---

## Stage 3 -- the distogram: training and leave-fold-out inference

**Training** (`core/predict.py:345 train_fold`): for fold f, `entries = peptides with fold != f`
(line 362) `+ fold_fragments(f)` (line 364) -> `dataset` (191): `X` stacked, labels
`Y = np.digitize(d, BIN_EDGES)` (line 200) from the true CA-CA distance of every pair. Measured
training sets: fold 0: 629 peptides + 6,001 fragments = 6,630 chains, 552,199 pairs; fold 1:
647 + 5,988 = 6,635 / 556,344; fold 2: 608 + 6,001 = 6,609 / 554,123; fold 3: 626 + 6,003 = 6,629 /
552,912; fold 4: 638 + 6,002 = 6,640 / 556,709. **The 787 peptides are 9.5% of the chains; the
6,003 shared fragments are the other 90.5%.**

Bins (`BIN_EDGES`, line 54): 16 edges 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0,
12.5, 14.0, 16.0, 19.0, 23.0 -> **17 bins**; `CENTRES` (line 57) = 4.0, 4.75, 5.25, 5.75, 6.25,
6.75, 7.25, 7.75, 8.5, 9.5, 10.5, 11.75, 13.25, 15.0, 17.5, 21.0, 25.0 A. Bin 0 is d < 4.5 (centre
4.0 by fiat), bin 16 is d >= 23 (centre 25.0 by fiat). Adjacent-centre gaps 0.5 .. 4.0 A.

Model (`MLP`, line 231): 183 -> 384 -> 384 -> 384 -> 17 with GELU, **372,881 parameters**
(measured), dropout 0 (`DROPOUT = 0.0`, line 335, with the reason at 315-334: dropout improves the
distance matrix and worsens the ranking), input standardisation (265-266). Loss (`fit`, 261):
soft-binned cross-entropy, target = Gaussian in BIN INDEX with sigma `smooth = 0.6` (line 269),
per-example weight 1 (`sample_w` is never passed by `train_fold`; `separation_weights`, line 208, is
unused on the deployed checkpoint), AdamW lr 2e-3, weight decay 1e-4, OneCycle, 40 epochs, batch
4096, torch seed 0. Checkpoint `distogram_models/fold{f}_esm_frag.pt` (`_model_path`, 339). One
model per fold (`N_ENSEMBLE = 1`, 336).

**Inference**: `core/pipeline.py:618 fold_model(fold)` -> `train_fold(fold, True, 5,
fragments=True)` loads the checkpoint (never trains; `guard_esm` at 535 refuses a training path);
`core/predict.py:426 Distogram.for_target(seq, model)` -> `prob = softmax(MLP(X))` **(npairs, 17)**;
`expected (npairs,)` = sum prob * CENTRES; `sd (npairs,)` (lines 406-408). Cached per target as
`s12/cache/disto_<pdb>.npz` (`s12/instrument.py:178 distogram`).

CREATED: the map from sequence to a distance distribution, learned from 6.6k chains. This is the
system's only sequence-conditioned structural information beyond retrieval. COMPRESSED: a
continuous distance -> 17 bins (log2 17 = 4.09 bits per pair at most; the L1 consumer in stage 4
takes only the median atom, so the effective per-pair target takes 17 values with a largest
adjacent gap of 4.0 A, S25 L6). Measured properties: over-confident by ~2x (z_sd 1.66-2.09, S25 L1),
24.1% of pairs multimodal, MAE ~2.0 A.

**COMMON-MODE ENTRY 1.** The predictor emits "a typical peptide of this length" (S18/S19): its
error is coherent across pairs and points the same way on every target (S19 L-series;
`error-shape-not-mae-decides-ranking`). Observable in principle only against a referent outside
the model: the native (ORACLE), or an independent structural source. It cannot be observed from the
posterior alone (calibrating its width or moving its location does not move the endpoint, S25 L2/L12).

---

## Stage 4 -- the risk table and the score

`core/predict.py:399 Distogram.__init__`, lines 409-423:

    w_p    = shell / (sd_p + 0.5)^g          shell = 1, g = 1  (score_weights.json absent; _score_weights:388)
    w_p   <- w_p / mean_p w_p                (line 417)
    grid   = 2.00, 2.05, ..., 39.95  A       (760 points, line 419)
    risk[p, t] = w_p * sum_c prob[p, c] * |grid_t - CENTRES_c|      (npairs, 760) float32, lines 420-422

`s12/instrument.py:200 shipped_score(dg, D)` (identical to `core/predict.py:443 score` and
`s7/debias.py:238 score_risk`): for a candidate's pair distances `D (npairs,)`, `g = clip(int((D -
2.0)/0.05), 0, 759)` (truncation toward zero), `score = mean_p risk[p, g_p]` -> one float, lower is
better. The distances are cast float32 -> float64 first (`s24/d_harness.py:259`), the reference
path's storage round trip.

TRANSFORMED: the posterior (17 numbers per pair) -> a 1-D per-pair risk curve whose minimiser is the
posterior MEDIAN (an L1 Bayes risk). COMPRESSED: (npairs,) candidate distances -> 1 scalar; the
0.05 A grid quantisation is harmless, the 17-atom quantisation of the target is not (S25 L6).
DESTROYED: which pairs disagree and in which direction (the signed deviation map; the set
transformer over it had a flat learning curve, S13), and every pair's contribution beyond its
weighted absolute deviation. The weight uses the width only through 1/(sd + 0.5).

---

## Stage 5 -- retrieval: the window universe and the K = 500 pool

Pinned universes: `s8/generate.py:159 stage_univ`, one `s8/generate_univ/<pdb>.npz` per target.
Library for a target of fold f: peptides with `fold != f and seq != target` + `fold_fragments(f)`
(lines 199-200); `_windows_all` (137) enumerates EVERY length-n window of every member:

    W     (nw, n, 3)  float32  A      CA coordinates of the window (the parent's real geometry)
    PHI   (nw, n)     float16  rad    the parent's torsions (float16: 10-bit mantissa, ~0.002 rad near pi)
    PSI   (nw, n)     float16  rad
    S     (nw, n)     int8            residue codes
    org   (nw,)       bool            True = peptide window, False = protein fragment window
    sim   (nw,)       float32         sum_i BLOSUM62[target_i, window_i]   (line 204)
    order (nw,)       int32           stable argsort of -sim                (line 205)
    rr    (nw,)       float32  A      ORACLE: CA-RMSD of each window to nat_ca (line 219)
    nat_ca (n, 3)     float32  A      ORACLE: model-1 native CA trace        (line 220)

Measured: nw = 7,016 .. 39,410 windows per target, median 17,088.

The pool: `s12/instrument.py:63 pool_idx` -> `order[:500]`; `s24/d_harness.py:185
Candidates.from_universe(pdb, k=500)` -> `W (500, n, 3)` float64 (from float32), `PHI/PSI (500, n)`
float64 (from float16), `nat_ca`, `oracle_rr` kept apart as labels, `meta.universe_idx`. Production
equivalent: `core/pipeline.py:672 retrieve` with `core/data.py:167 top_k` (STABLE argsort, ties
resolve by bank position; the unstable sort changes up to 47 of 500 members) and the float32 round
trip `_q` (666).

CREATED: **the only three-dimensional coordinates in the path** -- every candidate is a real
deposited backbone window. COMPRESSED: ~17k windows -> 500 by a sequence-only key that is an
integer BLOSUM sum over n positions with large tie sets; log2 C(17088, 500) = 3,252 bits of choice
resolved by a key whose Spearman with the true RMSD over the universe is +0.066
(`s8/generate.py` docstring, RESULT 2). DESTROYED: 97% of the universe; the parent's sequence
context beyond the window; the parent's identity except through `org`.

**COMMON-MODE ENTRY 2 (and the place it is measured).** S23 L9 (`s23/errdecomp.py`, n = 126): in
the medoid frame of the retained 75, mean_k |e_k|^2 = |ebar|^2 + mean_k |d_k|^2 holds exactly
(2.7e-14 A), with |ebar|^2 = 160.36 against mean |d_k|^2 = 63.82: **f = 0.676 of the pool's squared
error is a bias every member shares**, 50.7x the i.i.d. prediction 1/m, f > 1/m on 126/126. S24
L3/L5/L11: the shared bias is the SCORE's, not the corpus's and not BLOSUM's -- selection by the
prior (stages 4, 6, 7) installs it (L5's mechanism claim partly withdrawn in L11; the attribution
stands). Where it could in principle be observed: nowhere inside the pool (every within-pool
statistic is invariant to a common shift, by the identity above); only against a referent outside
the pool -- the native (ORACLE) or an independent native-free source that does not share the
prior's error. The one native-free contrast on record is the predicted-vs-realised Rg disagreement
(memory: `prediction-pool-disagreement-is-a-native-free-signal`; signal shown, Angstrom value not
measured).

---

## Stage 6 -- scoring the pool; the energy vector E

`s24/d_harness.py:250 score_shipped(cand)` -> `sc (500,)` via `I.pair_index`, `I.pair_dists`
(`s12/instrument.py:133`) and `shipped_score` (stage 4). `s25/phys_lib.py:144 channels(pdb)` asserts
`sc` bit-equal to `s24/cache_amber/<pdb>.npz :: score_dist` and the pool identity bit-equal to that
cache's `universe_idx` (lines 155-158). `s27/run_pool.py:99 channels_for` stores every channel as a
(500,) vector in `s27/cache/<pdb>.npz`; `DIS` is `sc`.

`E = zrank(sc)`: `s27/run_pool.py:68 zr` -> `s27/ham_lib.py:464 zrank` (= `s24/d_harness.py:274`, =
`core/pipeline.py:788 _zrank`): `(rankdata(x) - mean) / sd`. Output `E (500,)`, mean 0, sd 1.

COMPRESSED: 500 real scores -> 500 standardised ranks. The per-target information that survives is
the PERMUTATION (which candidate holds which rank), at most log2 500! = 3,767 bits; the SPECTRUM
{E_i} is the same ladder for every target up to tie-averaging (S25 L17: worst deviation 1.18% of
range across 8 targets, 0 of 8 exactly equal). In the production `quantum_stage`
(`core/pipeline.py:836-838`) E is formed over the top-128 in ascending order, so it is literally the
standardised ranks 1..128 for every target. DESTROYED: every score magnitude and every gap -- the
only thing left for a Hamiltonian to be "expressive about" is the label.

---

## Stage 7 -- the tie-safe top-m (the spine) and the candidate-basis encoding (the VQE arm)

Spine (`s27/run_vqe_chain.py:105-109, 125`): `key = rng_for(pdb, "tiekey").random(500)`
(`s27/run_pool.py:54`: sha256 of `"s27|<pdb>|tiekey"`), `order = np.lexsort((key, E))`, `top =
order[:75]` (`M = 75`, `s27/run_pool.py:38`). Exact ties (target 5H1H has an exact top-2 tie,
`s24/d_harness.py:419`) break by the random key, never by array order.

Output: a set of 75 indices. Information: log2 C(500, 75) = 300.6 bits of set choice out of the
3,767 bits of ordering; DESTROYED: the order inside the set and outside it, and every score.

VQE arm encoding (`s22/qcand_lib.py:120 Encoding(E)`): `n_qubits = ceil(log2 500) = 9`, `dim =
512`; `pad_energy = max(E) + 10 * sd(E)` (line 135) on the 12 padding states; `label = arange(512)`
(136; the gauge control `random_label`, 151, permutes it); `E (512,)` with the 500 real energies at
slots 0..499 (141-142). Production (`core/pipeline.py:836-838`): `n = 7`, the top-128 prefix, no
padding.

---

## Stage 8 -- the CVaR-VQE (research arm; NOT on the anchor's path)

`core/quantum.py:818 StatevectorCircuit(n, layers=3, ring=True)`:

    |psi(theta)> = prod_{L=1..3} [ U_ent . (x)_q RY(theta_{L,q}) ] |0..0>
    U_ent        = CNOT(0,1) CNOT(1,2) ... CNOT(n-2,n-1) CNOT(n-1,0), composed once into a basis permutation (`_entangler_permutation`, 841-856)
    theta        (3n,) rad :  27 at n = 9 (harness), 21 at n = 7 (production)
    psi          (2^n,) REAL amplitudes (RY and CNOT are real; `states_batch`, 859-882)
    p_theta      = psi^2 / sum   (`probs_batch`, 884)

`core/quantum.py:1023 run_cvar_vqe(E, alpha, T, n, layers, iters, restarts=1, seed, lr=0.15)`:
`theta0 ~ N(0, 0.6^2)` (1039), Adam beta 0.9/0.999, eps 1e-8 (1044-1046), `iters = 80` in the
harness (`s27/run_vqe_chain.py:39`, `s24/d_harness.py:321`), 50 in production (`Config.vqe_iters`,
`core/pipeline.py:183`). Objective `free_energy` (993):

    F(theta) = CVaR_alpha(E; p_theta) - T * H(p_theta),     H = -sum p log p   (clip 1e-15, 1014)
    CVaR_alpha = (1/alpha) [ sum_{E(x)<q} p(x) E(x) + (alpha - P(E<q)) q ]      (cvar_exact, 906-926)
    dF/dtheta_k = sum_x [dCVaR/dp(x) - T dH/dp(x)] . [p(x; theta_k + pi/2) - p(x; theta_k - pi/2)] / 2   (1018-1019)

Settings: harness alpha = 0.18, T = 0.5 (`s27/run_vqe_chain.py:39`; `s25/phys_lib.py:71-72`),
seeds 0 and 1; production `VQE_LFO` (`core/pipeline.py:113`): T = 0.3 on every fold, alpha = 1.0 on
folds 0, 3, 4 and 0.25 on folds 1, 2. Output: `p (2^n,)`, `cvar`, `H`.

CREATED: nothing from data. The inputs are E -- a target-independent ladder up to ties -- the
seed, alpha and T; p_theta is a function of those alone, so up to ties and the seed **the trained
distribution is the same distribution on every target** (S25 L17: "two trained states in the whole
deployment"). What differs per target is which candidate each basis index names. TRANSFORMED: a
27-parameter family on the 511-dimensional real sphere (dim so(512) = 130,816; the circuit has 27
generators) is optimised on a diagonal H whose exact simplex optimum at alpha = 1 is the Gibbs
state exp(-E/T)/Z; at the deployed settings the trained state sits 0.90 nats / 45% total variation
from it (S25 L15, `s25/QUANTUM.md` section 6.3). Cost per gradient: 2P = 54 (or 42) full
statevectors, batched.

---

## Stage 9 -- the tail readout

`s22/qcand_lib.py:299 exact_face(E, p, alpha)` -> `core/quantum.py:290 cvar_from_probs`: `order =
argsort(E, stable)` (306), `cum = cumsum(p[order])` (307), `take = clip(alpha - (cum - p[order]), 0,
p[order])` (308), `mass[order] = take`; `tail = mass > 0`, `tail_size`, `ess = alpha^2 / sum mass^2`.
`s22/qcand_lib.py:352 tail_candidates` maps tail bit-indices to candidate indices and drops
padding. `s24/d_harness.py:321 arm_vqe` returns `m = len(cands)`, `cands` sorted, `entropy_bits`,
`ess`. `gate_set_equality` (402) asserts the theorem: `cands` is a subset of the energy prefix that
reaches its own worst member (value-based, tie-safe).

DESTROYED, in order: (i) the **signs** of the 2^n real amplitudes (p = psi^2; 511 sign bits gone
before the objective is even formed); (ii) **the probabilities** -- once `tail` is a boolean mask the
mass vector is discarded and the next stage weights every survivor equally (uniform readout,
`readout_uniform`; the p-weighted readouts R2/R3 are +0.24 to +0.33 A WORSE, S28-L21); (iii) **the
order** inside the tail. What survives to the structure: the integer m (which rung of the classical
ladder the alpha-mass reaches) and, at exact zeros only, which prefix members are deleted. By the
set-equality theorem the survivors ARE `argsort(E)[:m]` (0 subset violations in 2,592 + 3,888 +
17,574 cells; `s25/QUANTUM.md` section 5), so the whole quantum stage reduces, for the structure,
to choosing m: 74 (seed 0) or 71 (seed 1) against the spine's 75 (S28-L21 anchors, `s27/results/vqe_rows.jsonl :: DIS`).

---

## Stage 10 -- the uniform coordinate average (the point cloud, 3.0483 A)

`s24/d_harness.py:288 readout_uniform(cand, idx)` -> `s12/instrument.py:119 coordinate_average(W[idx])`:

    P     = pairwise Kabsch CA-RMSD of the 75 members, (75, 75)     (pairwise_rmsd, 107; kabsch_rmsd_batch, 79)
    b     = argmin_k mean_j P[k, j]           the MEDOID              (medoid, 115)
    W'_k  = R_k W_k + t_k  superposed on W_b, reflection forbidden    (superpose_batch, 95-104)
    C     = mean_k W'_k                        (n, 3) A               (line 125)

Production equivalent `core/pipeline.py:924 average` (identical operator: `s8/consensus2.py:293
_medoid`, `superpose_batch`, `.mean(0)`).

Input 75 x n x 3 -> output n x 3: **compression 75x**. DESTROYED: the members' identities, every
pairwise relation among them (P is used only to pick the frame), the spread and covariance of the
ensemble, its multimodality (2-3 populated clusters per target, `core/pipeline.py:810`), and 22.3%
of the chain's length (mean virtual CA-CA 2.961 A against a native 3.812 A, worst bond 0.649 A,
`ARCHITECTURE.md` section 0). PASSED THROUGH UNCHANGED: the 68% common-mode error -- averaging
cancels only the idiosyncratic 32% (S23 L9: averaging is worth -0.655 A, 126/126, and cannot touch
|ebar|). Removing the four most deviant members costs +0.142 A while removing four at random costs
nothing (S23 L5): the diversity IS the variance reduction. Dimension: n x 3 - 6 rigid dof
(21..42). Point-cloud RMSD 3.0483 (`chain_rows.jsonl :: DIS :: rmsd_cloud`; not a structure).

---

## Stage 11 -- the ideal-geometry projection (the built chain, 3.2126 A)

`s24/d_harness.py:304 readout_projected` -> `s12/instrument.py:139 project(C, seq, fold, lam=0.3,
multi=True, maxiter=300)` -> `core/project.py:759 make_penalty("ramah")` and `933 lam_path(C, pen,
(0.0, 0.3), maxiter=300, multi=True, grad="exact")`:

    rung lam = 0   : fit_multi (912) -- L-BFGS-B from the FOUR generic starts (STARTS, 814: extended
                     (-120,130), alpha (-57,-47), beta (-139,135), PPII (-75,145) deg, every residue)
                     on  CA-RMSD( build(phi, psi), C ),  keep the lowest objective
    rung lam = 0.3 : fit_prior (890) warm-started from the lam = 0 solution on
                     CA-RMSD( build(phi, psi), C ) + 0.3 * pen(phi, psi),
                     PLUS fit_multi from the four starts again (multi=True, 955-958); keep the lower
    objective mode "exact" (_make_fg, 824-848): bit-exact builder build_ca_exact (374) and a ONE-SIDED
                     finite-difference gradient at FD_EPS = 1e-5 (185), 2n + 1 builds per evaluation

Penalty `RamaHingePenalty` (704): `mean_i max(0, t_class(i) - log P_fold(phi_i, psi_i))` over
interior residues, with P from the fold-disciplined count tables `s8/project_prior.json` (36 x 36
grid, wrapped-Gaussian smoothing sigma 1.5 bins, pseudo-count 0.5; `RB, SIGMA_BINS, PSEUDO`, 531-533)
and the hinge at the 5th percentile of real residues of that class (`HINGE_PCT`, 534; classes
gen/gly/pro/prepro, 536).

Geometry (149-155): N-CA 1.458, CA-C 1.525, C-N 1.329 A; N-CA-C 111.0, CA-C-N 116.2, C-N-CA 121.7
deg; omega = 180 deg -> virtual CA-CA 3.804 A, constant. Output: `phi, psi (n,)` rad and the CA trace
`(n, 3)` built from them (`_emit`, 817, via the reference builder), plus `fit_ca` (the lam = 0 arm).

Input n x 3 (3n - 6 dof) -> 2n torsions (the first phi and last psi do not move a CA; effectively
2n - 2) -> n x 3 on the ideal manifold. DESTROYED: **the cloud's distances** -- every virtual bond is
reset to 3.804 A, so the contraction is undone by construction and cis-peptides (~2.9 A) cannot be
represented; the cloud's exact shape (the nearest manifold point is taken). The problem is DEGENERATE:
two torsion branches sit at near-equal objective and a 1e-13 A change in the forward map routes some
L-BFGS-B trajectories to the other branch (`core/project.py` module docstring: the reference disagrees
with itself under a rigid motion by up to 1.6 A on some targets), which is why the shipped mode is
bit-exact rather than merely accurate. Cost: +0.1664 A over the cloud (median +0.098; the chain is
better on 16 of 126; rho with the target's own contraction +0.518; `ARCHITECTURE.md` section 0).
Built-chain RMSD 3.2126 (`chain_rows.jsonl :: DIS :: rmsd_chain`, n = 126, reproduced 2026-09-19).

---

## Stage 12 -- optional restrained AMBER relaxation (NOT on the anchor's path)

`core/pipeline.py:989 relax` -> `core/amber.py:1497 refine_coords(seq, rep, coords, k_restraint=10,
steps=0)`: ff14SB/GBn2 via OpenMM, positional restraints on N, CA, C (`BACKBONE_ATOMS`, 1260) at
k = 10 kcal/mol/A^2 (`Config.amber_k`, `core/pipeline.py:129`), minimised to convergence (`steps = 0`,
130), convergence gate 1000 kcal/mol (`CONVERGE_MAX_KCAL`, 1282). Production record `rmsd_full`
3.236 mean (`STATE_BRIEF` section 2), i.e. +0.02 A over the built chain; a validity stage. The
chain_rows basis does not apply it.

---

## Stage 13 -- evaluation

`s12/instrument.py:91 ca_rmsd(a, b)` = `kabsch_rmsd_batch` (79-88): centre both, SVD of the
cross-covariance, reflection corrected by the sign of det, RMSD = sqrt((|A|^2 + |B|^2 - 2 sum S)/n).
Native = the universe's `nat_ca` (stage 5; model 1 of the deposited PDB, float32 stored, float64
read). Rows: `s27/run_vqe_chain.py:128-130` -> `rmsd_cloud = ca_rmsd(C, nat_ca)`, `rmsd_chain =
ca_rmsd(ca, nat_ca)` per (config, pdb). Statistics: `s24/stats_lib.py:79 compare(a, b, folds)`:
d = a - b per target, SE = sd(d)/sqrt(n), MDE = 2.8016 x SE (`MDE_K`, 56), 4,000-draw iid
bootstrap and a cluster bootstrap that resamples the 5 folds with replacement (124-127), W/L/ties,
retrodesign power / Type-M / Type-S, verdict (`_verdict`, 145: |effect| must exceed its own MDE AND
the fold CI must exclude zero). Folds: `pinned_folds` (65) reads the universes' `fold`.

---

## Summary table: what is created, kept and lost

| stage | object in -> out | dimension / entropy in -> out | created | destroyed |
|---|---|---|---|---|
| 0 corpus | deposits -> 787 + 6,003 chains | -- | sequence-structure pairs | side chains, context, other NMR models |
| 1 ESM | seq (n log2 20 = 39..69 bits) -> (n,1280)+(n,n) -> (n,32) | 1280 -> 32 per residue (40x) | language-model context | 1248 of 1280 directions |
| 2 features | -> (npairs, 183) | 183 per pair | -- | j-i = 1 pairs; context beyond +-2 |
| 3 distogram | (npairs,183) -> (npairs,17) | <= 4.09 bits per pair | seq -> distance law (6.6k chains) | continuous distance (17 atoms, gaps to 4 A) |
| 4 score | (npairs,17) -> risk (npairs,760); D (npairs,) -> 1 scalar | npairs -> 1 | -- | signed deviation map; width beyond 1/(sd+0.5) |
| 5 retrieval | ~17k windows -> 500 | log2 C(17088,500) = 3,252 bits of choice by a key at rho +0.066 | the ONLY 3-D coordinates | 97% of the universe; parent context |
| 6 energy | 500 scores -> 500 ranks | 3,767 bits of permutation; spectrum target-independent | -- | score magnitudes and gaps |
| 7 top-75 / encoding | 500 -> set of 75 (or 512-state register) | 300.6 bits of set choice | -- | order; scores |
| 8 CVaR-VQE | E (512,) -> p_theta (512,) | 27 params on a 511-sphere; same E on every target | nothing from data | -- |
| 9 tail readout | p -> mask -> m candidates | 512 reals -> 1 integer m | -- | amplitude signs (511 bits), probabilities, order |
| 10 average | 75 x n x 3 -> n x 3 | 75x compression; 3n-6 dof | -- | identities, pairwise relations, spread, modes, 22% of length; common-mode error passed through |
| 11 projection | n x 3 -> 2n torsions -> n x 3 | 3n-6 -> 2n-2 | ideal covalent geometry | the cloud's distances; branch chosen by arithmetic |
| 12 AMBER | n x 3 -> n x 3 (optional) | -- | all-atom validity | ~0.02 A |
| 13 RMSD | n x 3 vs nat_ca -> 1 number | 3n-6 -> 1 | -- | everything but the Kabsch residual |

## Where the 68% common-mode error enters and where it could be seen

- **Enters at stage 3** as the predictor's coherent "typical peptide" error and **at stage 5** as
  the pool members' shared bias; **stages 4/6/7 transmit it** (selection by the prior installs the
  prior's error in the retained set: S24 L5, attribution retained after L11; confidence weighting
  reduces the transmitted share by ~6% and buys ~0 A: S24 L6); **stage 10 cannot remove it**
  (exact identity, S23 L9).
- **Measured at**: stage 10's output, ORACLE, S23 L9 (`s23/results/errdecomp.json`): f = 0.676.
- **Could in principle be observed, native-free**: only by a contrast between two objects that do
  not share the error -- (a) the prior's prediction against the pool's realised geometry (predicted
  vs realised Rg, memory `prediction-pool-disagreement-is-a-native-free-signal`); (b) the
  sequence-conditioned answer against a sequence-blind answer (S29 H1, lane O's probe; the record's
  geometry, S24 L2/L3, expects the sign to be wrong); (c) any external structural source. It cannot
  be seen from any within-pool statistic, from the posterior's width, or from the VQE's state
  (which sees only the ladder).

## Which stages the anchor and the arms actually traverse

    chain_rows.jsonl :: DIS (3.2126)   : 0 1 2 3 4 5 6 7(top-75) 10 11 13
    chain_rows.jsonl :: rmsd_cloud     : 0 1 2 3 4 5 6 7(top-75) 10 13
    vqe_rows.jsonl :: DIS (S27/S28)    : 0 1 2 3 4 5 6 7(encoding) 8 9 10 13       (point cloud)
    production core.pipeline record    : 0 1 2 3 4 5 6 7(top-75) 10 11 12 13; quantum: null
    S25 seven-configuration suite      : as vqe_rows, with E = zrank(sum zrank(channels))
