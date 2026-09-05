"""Inverse folding as a selection objective: score the SEQUENCE given the STRUCTURE.

THE QUESTION THIS ANSWERS
-------------------------
Sprint 7 finding 5 measured that the shipped distance objective is MISSPECIFIED: on a
126-target instrument the native is the objective's argmin on 3/126 targets and sits at
the 36.8th percentile of its own pool's score distribution. Finding 10 measured the same
for all-atom Amber and got a WORSE answer: 0/70 argmin, 54.3rd mean percentile (median
61.5), and significantly worse than the distance prior once radius of gyration is
partialled out.

Both of those channels run the arrow the same way -- predict something about the structure
from the sequence, then measure agreement. This module runs it the other way.

    score(C) = -(1/n) * sum_i log P( s_i | local backbone geometry of C at i )

where `s` is the TARGET'S ACTUAL SEQUENCE and `C` is a candidate structure. The native is
the structure the real sequence actually folds into, so it should explain that sequence
better than a decoy does. Lower is better, like every other term in this project.

WHY IT MIGHT WORK WHERE THE OTHER TWO FAILED
  * Inverse folding is far better posed than forward folding. Many sequences fold to one
    structure, so structure -> residue propensity is a much tighter mapping than
    sequence -> structure. ProteinMPNN reaches ~52% native sequence recovery; nothing in
    this project predicts structure from sequence at anything like the corresponding
    quality.
  * It should survive the distribution shift that killed finding 2 (30x more structural
    data made the distance predictor MONOTONICALLY worse, because a fragment inside a
    folded protein has a more compact conformational distribution than an isolated
    peptide). Local geometry -> residue propensity is close to universal; isolated-peptide
    distance distributions are not.
  * It is not compactness-driven, so it cannot fail Amber's way. (Verified anyway: every
    reported correlation is also reported with radius of gyration partialled out.)
  * It is CHIRALITY-SENSITIVE. Finding 4 proved reflection leaves the CA-CA distance
    matrix EXACTLY unchanged. The features here are signed volumes in a local frame, so a
    mirrored structure scores differently. `test_invfold.py` asserts this.

GEOMETRY PROVENANCE -- the trap finding 10 documents
----------------------------------------------------
Pool candidates are ideal-geometry backbones rebuilt from a real fragment's (phi, psi) via
`geo.build_backbone_batch`: fixed bond lengths, fixed bond angles, fixed trans omega. A
deposited native is an experimental structure. Scoring one against the other compares
geometry PROVENANCE, not conformation.

So this module holds provenance constant in BOTH directions:
  * TRAINING examples are length-9..16 windows of `prots/` chains, REBUILT from the
    window's own (phi, psi) through the same builder. The model therefore sees exactly the
    kind of object it will be asked to score -- an isolated ideal-geometry short backbone
    -- and never a deposited one.
  * The NATIVE, when scored, is rebuilt from its own (phi, psi) through the same builder.

FEATURES -- geometry only, NEVER the residue's own identity
-----------------------------------------------------------
Per residue i, all rigid-motion invariant:
  * backbone torsions of i-1, i, i+1 as sin/cos (12)
  * for each of the k nearest CA neighbours j (by CA-CA distance):
      - 8 radial basis functions of d(i, j)
      - the unit vector from i to j expressed in i's own N/CA/C frame (3) -- SIGNED, this
        is where the chirality sensitivity lives
      - the relative ORIENTATION of j's frame in i's, R_i^T R_j, as its first two rows
        (the continuous 6D rotation representation) (6).
        A distance matrix throws this away entirely; it is the main thing this
        representation has that the distance channel structurally cannot have.
      - a 12-way one-hot of bucketed sequence separation, plus the signed separation (13)
  * neighbour-count densities at 6/8/10 A, normalised radial position, position in chain,
    chain length (7)

Variant `nbid` additionally supplies each neighbour's amino-acid identity as a 20-way
one-hot. That is strictly more information and is legitimate at scoring time (we know the
target sequence; we are asking which STRUCTURE explains it). But it risks the model
learning sequence-composition statistics rather than geometry, so both variants are built
and both are reported. The residue's OWN identity is never an input in either variant.

LEAKAGE
-------
Training reads `prots/*.pdb` only. `prots/` has zero PDB-id overlap with the peptide DB,
the 24-target dev set, or the 60-target benchmark (verified by `stage_leak`, which also
re-verifies it here rather than trusting the record). On top of that the same identity
filter `distogram._fold_fragments` applies is used at TWO levels: whole chains above the
threshold to any peptide-DB sequence are dropped, and then every extracted WINDOW above
the threshold to any peptide-DB sequence is dropped. No native coordinate, distance,
torsion, contact or RMSD of any evaluation target enters training or scoring.

DIAGNOSTIC vs DEPLOYABLE
------------------------
  * DEPLOYABLE: `features`, `InvFoldNet`, `SelfConsistency.score`. These read a candidate
    structure and the target sequence. No native.
  * ORACLE DIAGNOSTIC: `stage_native` / the `native_pct`, `native_is_argmin`, `rr`,
    `rho*` quantities in `stage_select`. These deliberately read natives to answer whether
    the objective is correctly specified, and to report RMSD after a selection that was
    made without them. Nothing deployable imports them.

STAGES
------
    python -m s8.invfold leak      # leakage audit (run first, it gates everything)
    python -m s8.invfold scan      # parse prots/ -> torsion cache      (~8 min)
    python -m s8.invfold windows   # leak-filtered training windows     (~3 min)
    python -m s8.invfold train     # both variants                      (~25 min)
    python -m s8.invfold native    # THE DECISIVE EXPERIMENT
    python -m s8.invfold select    # paired selection vs the distogram
    python -m s8.invfold report
"""
import glob
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np

import peptide_db as pdb
import protein_geometry as geo

HERE = os.path.dirname(os.path.abspath(__file__))
SCAN_CACHE = os.path.join(HERE, "invfold_scan.npz")
WIN_CACHE = os.path.join(HERE, "invfold_windows.npz")
MODEL_DIR = os.path.join(HERE, "invfold_models")

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA)}
NAA = len(AA)

K_NEIGH = 8
N_RBF = 8
RBF_MU = np.linspace(3.0, 16.0, N_RBF).astype(np.float32)
RBF_SIG = float(RBF_MU[1] - RBF_MU[0])
#: Bucketed rather than raw: what matters is "adjacent", "one turn away", "far in
#: sequence but close in space", not the exact integer.
SEP_EDGES = np.array([-16, -8, -4, -2, -1, 0, 1, 2, 4, 8, 16])
N_SEP = len(SEP_EDGES) + 1          # 12

MIN_LEN, MAX_LEN = 9, 16
PER_NB = N_RBF + 3 + 6 + N_SEP + 1          # 30
N_TOR, N_GLOB = 12, 7


def n_features(nbid: bool = False) -> int:
    per = PER_NB + (NAA if nbid else 0)
    return K_NEIGH * per + N_TOR + N_GLOB


# ===================================================================== features
def _frames(N, CA, C):
    """Per-residue right-handed orthonormal frame ``(B, n, 3, 3)``, rows = (e1, e2, e3).

    e1 along CA->C, e2 the N-direction orthogonalised against it, e3 = e1 x e2. The cross
    product is what makes the frame -- and therefore every neighbour coordinate expressed
    in it -- change sign under reflection.
    """
    e1 = C - CA
    e1 = e1 / np.maximum(np.linalg.norm(e1, axis=-1, keepdims=True), 1e-9)
    v = N - CA
    e2 = v - (v * e1).sum(-1, keepdims=True) * e1
    e2 = e2 / np.maximum(np.linalg.norm(e2, axis=-1, keepdims=True), 1e-9)
    e3 = np.cross(e1, e2)
    return np.stack([e1, e2, e3], axis=-2)


def _knn(d, n, k):
    """Indices of the k nearest CA neighbours, ties broken DETERMINISTICALLY by index.

    Not a detail. Candidates are ideal-geometry rebuilds, so consecutive CA-CA distances
    are all EXACTLY 3.80 A and a residue's two sequence neighbours tie exactly. Under a
    rigid rotation those distances move by ~1e-15, so an unstable sort silently returns a
    DIFFERENT neighbour set and the score stops being rotation-invariant --
    `test_score_invariant_to_rigid_motion` catches it. Rounding to 1e-6 A (far below any
    physical resolution) and using a stable sort makes the neighbour graph exact.
    """
    key = np.round(d + np.eye(n)[None] * 1e6, 6)
    return np.argsort(key, axis=-1, kind="stable")[:, :, :k]


def features(coords: Dict[str, np.ndarray], phi=None, psi=None,
             nb_ids: Optional[np.ndarray] = None) -> np.ndarray:
    """``(B, n, D)`` rigid-motion-invariant per-residue geometry features.

    DEPLOYABLE. Reads only the candidate structure (and, for the `nbid` variant, the
    target sequence which is a given of the problem). The residue's own identity is never
    included -- that is the label.

    `coords` holds N, CA, C as ``(B, n, 3)`` or ``(n, 3)``. `nb_ids` is ``(n,)`` integer
    amino-acid indices; pass it only for the `nbid` variant.
    """
    N, CA, C = coords["N"], coords["CA"], coords["C"]
    if CA.ndim == 2:
        N, CA, C = N[None], CA[None], C[None]
    N = np.asarray(N, np.float64); CA = np.asarray(CA, np.float64); C = np.asarray(C, np.float64)
    B, n, _ = CA.shape
    F = _frames(N, CA, C)                                    # (B, n, 3, 3)

    diff = CA[:, None, :, :] - CA[:, :, None, :]             # (B, i, j, 3): j - i
    d = np.linalg.norm(diff, axis=-1)
    k = min(K_NEIGH, max(1, n - 1))
    idx = _knn(d, n, k)

    bi = np.arange(B)[:, None, None]
    ni = np.arange(n)[None, :, None]
    dk = np.take_along_axis(d, idx, axis=-1)                 # (B, n, k)
    vk = diff[bi, ni, idx]                                   # (B, n, k, 3)

    rbf = np.exp(-((dk[..., None] - RBF_MU) / RBF_SIG) ** 2).astype(np.float32)
    # Neighbour direction in residue i's OWN frame. Signed -> chirality-sensitive.
    local = np.einsum("bnij,bnkj->bnki", F, vk) / np.maximum(dk[..., None], 1e-9)
    # Relative ORIENTATION of j's frame in i's, R_i^T R_j -- the SE(3) information a
    # distance matrix cannot carry at all. Kept as the first two ROWS of the matrix (the
    # standard continuous 6D rotation representation) rather than a quaternion: the
    # quaternion sign convention w >= 0 is DISCONTINUOUS at w = 0, which breaks exact
    # rigid-motion invariance of the features. `test_score_invariant_to_rigid_motion`
    # catches that, and it is not a test artefact -- it would be a real representational
    # discontinuity in the deployed score.
    # F holds basis vectors as ROWS, so R_i = F_i^T and R_i^T R_j = F_i F_j^T.
    Fj = F[bi, idx]                                          # (B, n, k, 3, 3)
    Rrel = np.einsum("bnax,bnkcx->bnkac", F, Fj)
    q = Rrel[..., :2, :].reshape(Rrel.shape[:-2] + (6,)).astype(np.float32)

    sep = (idx - np.arange(n)[None, :, None]).astype(np.float32)
    sb = np.digitize(sep, SEP_EDGES)
    sep_oh = np.zeros(sb.shape + (N_SEP,), np.float32)
    np.put_along_axis(sep_oh, sb[..., None], 1.0, axis=-1)

    parts = [rbf, local.astype(np.float32), q, sep_oh, np.tanh(sep / 8.0)[..., None]]
    if nb_ids is not None:
        oh = np.zeros((NAA + 1,), np.float32)
        table = np.eye(NAA, dtype=np.float32)
        lab = np.asarray(nb_ids, int)
        parts.append(table[lab][idx])                        # (B, n, k, 20)
        del oh
    nb = np.concatenate(parts, axis=-1)
    if k < K_NEIGH:                                          # pad short chains
        nb = np.concatenate(
            [nb, np.zeros((B, n, K_NEIGH - k, nb.shape[-1]), np.float32)], axis=2)
    nb = nb.reshape(B, n, -1)

    if phi is None or psi is None:
        tor = np.zeros((B, n, N_TOR), np.float32)
    else:
        ph = np.atleast_2d(np.asarray(phi, float))
        ps = np.atleast_2d(np.asarray(psi, float))
        cols = []
        for off in (-1, 0, 1):
            j = np.clip(np.arange(n) + off, 0, n - 1)
            cols += [np.cos(ph[:, j]), np.sin(ph[:, j]),
                     np.cos(ps[:, j]), np.sin(ps[:, j])]
        tor = np.stack(cols, axis=-1).astype(np.float32)
        if tor.shape[0] == 1 and B > 1:
            tor = np.broadcast_to(tor, (B, n, N_TOR))

    rad = np.linalg.norm(CA - CA.mean(1, keepdims=True), axis=-1)
    ar = np.arange(n)
    glob = np.stack([
        (d < 6.0).sum(-1) - 1.0, (d < 8.0).sum(-1) - 1.0, (d < 10.0).sum(-1) - 1.0,
        rad / n ** 0.38,
        np.broadcast_to(ar / max(1, n - 1), (B, n)),
        np.broadcast_to(np.minimum(ar, n - 1 - ar) / max(1, n), (B, n)),
        np.full((B, n), n / 26.0),
    ], axis=-1).astype(np.float32)

    return np.concatenate([nb, tor, glob], axis=-1).astype(np.float32)


# ===================================================================== leakage
def held_sequences() -> List[str]:
    """Every sequence that must never be learnable: the whole peptide DB.

    Deliberately broader than "the evaluation targets". The 126-target tuning instrument,
    the 24-target dev set and the 60-target benchmark are all drawn from `peptide_db`, so
    excluding the entire DB makes ONE model safe for every target and removes any
    possibility of a fold-bookkeeping mistake. It is strictly stricter than the per-fold
    discipline `distogram._fold_fragments` applies.
    """
    return sorted({p.seq for p in pdb.load()})


class _KmerIndex:
    """3-mer inverted index over held sequences: the same prefilter `_fold_fragments` uses,
    made O(shared kmers) instead of O(n_held) per query."""

    def __init__(self, seqs: Sequence[str]):
        self.seqs = list(seqs)
        self.index: Dict[str, List[int]] = {}
        for i, s in enumerate(self.seqs):
            for km in pdb._kmers(s):
                self.index.setdefault(km, []).append(i)

    def leaks(self, seq: str, threshold: float = pdb.IDENTITY_THRESHOLD) -> bool:
        cand = set()
        for km in pdb._kmers(seq):
            cand.update(self.index.get(km, ()))
        for i in cand:
            if pdb.identity(self.seqs[i], seq) >= threshold:
                return True
        return False


def stage_leak():
    """Audit, run before anything else. Writes `s8/invfold_leak.json`."""
    prot_ids = sorted({os.path.basename(p)[:-4].upper() for p in glob.glob("prots/*.pdb")})
    pep_ids = {p.pdb.upper() for p in pdb.load()}
    dev_ids = {p.pdb.upper() for p in pdb.dev_set(24)}
    bench_ids = {p.pdb.upper() for p in pdb.benchmark()}
    out = {
        "n_prots_files": len(prot_ids),
        "n_peptide_db": len(pep_ids),
        "overlap_prots_peptidedb": sorted(set(prot_ids) & pep_ids),
        "overlap_prots_dev24": sorted(set(prot_ids) & dev_ids),
        "overlap_prots_bench60": sorted(set(prot_ids) & bench_ids),
        "identity_threshold": pdb.IDENTITY_THRESHOLD,
        "n_held_sequences": len(held_sequences()),
    }
    for k in ("overlap_prots_peptidedb", "overlap_prots_dev24", "overlap_prots_bench60"):
        out["n_" + k] = len(out[k])
    json.dump(out, open(os.path.join(HERE, "invfold_leak.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if not isinstance(v, list)}, indent=1))
    return out


# ===================================================================== data
def stage_scan(limit: Optional[int] = None):
    """Parse every `prots/*.pdb` to (sequence, phi, psi). Cached, ~8 min."""
    if os.path.exists(SCAN_CACHE):
        print(f"{SCAN_CACHE} exists")
        return
    paths = sorted(glob.glob("prots/*.pdb"))
    if limit:
        paths = paths[:limit]
    seqs, phis, psis, ids = [], [], [], []
    t0 = time.time()
    for k, p in enumerate(paths):
        try:
            seq, coords, phi, psi = geo.native_coords_from_pdb(p)
        except Exception:                                            # noqa: BLE001
            continue
        phi = np.asarray(phi, np.float32); psi = np.asarray(psi, np.float32)
        if len(seq) < MIN_LEN or len(phi) != len(seq) or len(psi) != len(seq):
            continue
        if not (np.isfinite(phi).all() and np.isfinite(psi).all()):
            continue
        seqs.append(seq); phis.append(phi); psis.append(psi)
        ids.append(os.path.basename(p)[:-4].upper())
        if (k + 1) % 2000 == 0:
            print(f"  {k+1}/{len(paths)} kept {len(seqs)} "
                  f"({time.time()-t0:.0f}s)", flush=True)
    np.savez_compressed(SCAN_CACHE + ".tmp.npz",
                        seqs=np.array(seqs, dtype=object),
                        phis=np.array(phis, dtype=object),
                        psis=np.array(psis, dtype=object),
                        ids=np.array(ids, dtype=object))
    os.replace(SCAN_CACHE + ".tmp.npz", SCAN_CACHE)
    print(f"scanned {len(paths)} files, kept {len(seqs)} chains, "
          f"{sum(len(s) for s in seqs)} residues in {time.time()-t0:.0f}s")


def stage_windows(per_chain: int = 30, seed: int = 0):
    """Leak-filtered length-9..16 training windows. Cached.

    The filter is applied at WINDOW level with the `distogram._fold_fragments` rule
    (identity >= 0.6 to any held-out sequence, with the same 3-mer prefilter). Window
    level is the right unit and is exactly the production discipline: `_fold_fragments`
    filters FRAGMENTS, which are windows, and a window is what this model trains on. A
    whole-chain filter was measured first and abandoned -- against a 9..16mer a
    170-residue chain needs a gapped alignment per candidate, ~1.1 ms each, and at 13,751
    chains x ~200 kmer-matched candidates that is over an hour to remove sequences the
    window filter removes anyway.
    """
    per_chain, seed = int(per_chain), int(seed)
    if os.path.exists(WIN_CACHE):
        print(f"{WIN_CACHE} exists")
        return
    z = np.load(SCAN_CACHE, allow_pickle=True)
    seqs, phis, psis = list(z["seqs"]), list(z["phis"]), list(z["psis"])
    ki = _KmerIndex(held_sequences())
    rng = np.random.default_rng(seed)
    keep_chain = list(range(len(seqs)))

    W_seq, W_phi, W_psi, W_len = [], [], [], []
    n_win_drop = 0
    t0 = time.time()
    for c, i in enumerate(keep_chain):
        s, ph, ps = seqs[i], phis[i], psis[i]
        m = len(s)
        if m < MIN_LEN:
            continue
        for _ in range(per_chain):
            L = int(rng.integers(MIN_LEN, min(MAX_LEN, m) + 1))
            if m - L < 0:
                continue
            a = int(rng.integers(0, m - L + 1))
            w = s[a:a + L]
            if any(ch not in AA_INDEX for ch in w):
                continue
            if ki.leaks(w):
                n_win_drop += 1
                continue
            W_seq.append(np.array([AA_INDEX[ch] for ch in w], np.int8))
            W_phi.append(ph[a:a + L]); W_psi.append(ps[a:a + L]); W_len.append(L)
        if (c + 1) % 2000 == 0:
            print(f"  {c+1}/{len(keep_chain)} chains, {len(W_seq)} windows "
                  f"({time.time()-t0:.0f}s)", flush=True)
    np.savez_compressed(WIN_CACHE + ".tmp.npz",
                        seq=np.array(W_seq, dtype=object),
                        phi=np.array(W_phi, dtype=object),
                        psi=np.array(W_psi, dtype=object),
                        L=np.array(W_len, np.int16),
                        n_chains_kept=len(keep_chain), n_chains_total=len(seqs),
                        n_window_leaks=n_win_drop, per_chain=per_chain)
    os.replace(WIN_CACHE + ".tmp.npz", WIN_CACHE)
    print(f"{len(W_seq)} windows ({sum(W_len)} residues), {n_win_drop} window-level "
          f"leaks dropped, in {time.time()-t0:.0f}s")


def build_xy(nbid: bool, max_res: int = 1_500_000, seed: int = 0, verbose: bool = True):
    """Feature matrix and labels from the cached windows.

    Windows are grouped by length so `build_backbone_batch` and `features` run batched.
    Every window is REBUILT from its own (phi, psi), so training geometry provenance is
    identical to candidate geometry provenance.
    """
    z = np.load(WIN_CACHE, allow_pickle=True)
    seq, phi, psi, L = z["seq"], z["phi"], z["psi"], z["L"].astype(int)
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(L))
    cum = np.cumsum(L[order])
    take = order[cum <= max_res] if cum[-1] > max_res else order
    # Preallocated, not accumulated-and-concatenated: at 1.5 M residues x 419 features the
    # concatenate would briefly hold two 2.5 GB copies, and this box has to fit two heavy
    # jobs at once.
    n_res = int(L[take].sum())
    D = n_features(nbid)
    X = np.empty((n_res, D), np.float32)
    Y = np.empty(n_res, np.int64)
    w = 0
    t0 = time.time()
    for ell in range(MIN_LEN, MAX_LEN + 1):
        sel = take[L[take] == ell]
        if len(sel) == 0:
            continue
        for s in range(0, len(sel), 4096):
            b = sel[s:s + 4096]
            ph = np.stack([phi[i] for i in b]).astype(float)
            ps = np.stack([psi[i] for i in b]).astype(float)
            sq = np.stack([seq[i] for i in b]).astype(int)
            # Same end-torsion convention as `build_pool`, so the model is trained on the
            # convention it will be scored with. A training window sliced from the middle
            # of a chain has real torsions at both ends; a scored peptide never does.
            ph[:, 0] = geo.DEFAULT_PHI
            ps[:, -1] = geo.DEFAULT_PSI
            c = geo.build_backbone_batch(ph, ps)
            f = _features_batch(c, ph, ps, sq if nbid else None)
            k = f.shape[0] * f.shape[1]
            X[w:w + k] = f.reshape(k, D)
            Y[w:w + k] = sq.reshape(k)
            w += k
        if verbose:
            print(f"  L={ell}: {len(sel)} windows, {w} residues "
                  f"({time.time()-t0:.0f}s)", flush=True)
    assert w == n_res, (w, n_res)
    return X, Y


def _features_batch(coords, phi, psi, seqs_int: Optional[np.ndarray]):
    """`features` for a batch whose neighbour identities differ per row.

    `features` takes a single `nb_ids` vector because at scoring time every candidate
    shares the target's sequence. In training each window has its own, so the one-hot
    block is assembled here instead.
    """
    f = features(coords, phi, psi, nb_ids=None)
    if seqs_int is None:
        return f
    B, n, _ = f.shape
    N, CA, C = coords["N"], coords["CA"], coords["C"]
    d = np.linalg.norm(CA[:, None, :, :] - CA[:, :, None, :], axis=-1)
    k = min(K_NEIGH, max(1, n - 1))
    idx = _knn(d, n, k)
    table = np.eye(NAA, dtype=np.float32)
    oh = table[seqs_int][np.arange(B)[:, None, None], idx]        # (B, n, k, 20)
    if k < K_NEIGH:
        oh = np.concatenate([oh, np.zeros((B, n, K_NEIGH - k, NAA), np.float32)], axis=2)
    # Splice the identity block into each neighbour slot so the layout matches
    # `features(..., nb_ids=...)` exactly.
    nb = f[..., :K_NEIGH * PER_NB].reshape(B, n, K_NEIGH, PER_NB)
    rest = f[..., K_NEIGH * PER_NB:]
    return np.concatenate(
        [np.concatenate([nb, oh], -1).reshape(B, n, -1), rest], -1).astype(np.float32)


# ===================================================================== model
class InvFoldNet:
    """MLP over local backbone geometry -> 20-way residue identity. DEPLOYABLE."""

    def __init__(self, d_in: int, width: int = 512, depth: int = 3,
                 dropout: float = 0.1, seed: int = 0):
        import torch
        import torch.nn as nn
        torch.set_num_threads(int(os.environ.get("NT", "2")))
        torch.manual_seed(seed)
        self.torch = torch
        self.d_in = d_in
        drop = [nn.Dropout(dropout)] if dropout > 0 else []
        layers = [nn.Linear(d_in, width), nn.GELU()] + list(drop)
        for _ in range(depth - 1):
            layers += [nn.Linear(width, width), nn.GELU()] + list(drop)
        layers += [nn.Linear(width, NAA)]
        self.net = nn.Sequential(*layers)
        self.width, self.depth = width, depth
        self.mu = np.zeros(d_in, np.float32)
        self.sd = np.ones(d_in, np.float32)
        self.bg = np.full(NAA, -math.log(NAA), np.float32)   # background log P(aa)

    def fit(self, X, Y, epochs: int = 12, batch: int = 8192, lr: float = 2e-3,
            val_frac: float = 0.05, verbose: bool = True):
        """MUTATES `X` (standardised in place) to keep peak memory at one copy.

        CAVEAT on `val_recovery`. The split is by RESIDUE, and one `prots/` chain
        contributes ~30 overlapping windows, so the same residue's geometry can appear on
        both sides. It is a TRAINING diagnostic and NOT a generalisation estimate -- it
        will read far too high. The honest held-out recovery is `native_recovery` in the
        selection pass, measured on the 126 tuning targets' own natives, which are
        excluded from training by sequence identity.
        """
        t = self.torch
        n = len(X)
        rs = np.random.default_rng(0).permutation(n)
        nv = int(n * val_frac)
        vi, ti = rs[:nv], rs[nv:]
        # Standardise IN PLACE and hand torch the same buffer. X is ~1.5 GB at the
        # production size; the obvious `t.as_tensor((X - mu) / sd)` would hold three
        # copies of it at once and this box has to fit two heavy jobs.
        sub = ti[:200_000]
        self.mu, self.sd = X[sub].mean(0), X[sub].std(0) + 1e-3
        cnt = np.bincount(Y[ti], minlength=NAA).astype(np.float64) + 1.0
        self.bg = np.log(cnt / cnt.sum()).astype(np.float32)
        X -= self.mu
        X /= self.sd
        Xn = t.from_numpy(X)
        Yt = t.as_tensor(np.ascontiguousarray(Y, np.int64))
        opt = t.optim.AdamW(self.net.parameters(), lr=lr, weight_decay=1e-4)
        per = int(np.ceil(len(ti) / batch))
        sched = t.optim.lr_scheduler.OneCycleLR(opt, max_lr=lr,
                                                total_steps=epochs * per + 1)
        lossf = t.nn.CrossEntropyLoss(label_smoothing=0.05)
        g = t.Generator().manual_seed(0)
        self.history = []
        for ep in range(epochs):
            self.net.train()
            perm = t.as_tensor(ti)[t.randperm(len(ti), generator=g)]
            tot = 0.0
            for s in range(0, len(ti), batch):
                b = perm[s:s + batch]
                loss = lossf(self.net(Xn[b]), Yt[b])
                opt.zero_grad(); loss.backward(); opt.step(); sched.step()
                tot += float(loss.detach()) * len(b)
            rec = {"epoch": ep, "train_loss": tot / max(len(ti), 1)}
            if nv:
                self.net.eval()
                with t.no_grad():
                    lo = t.log_softmax(self.net(Xn[t.as_tensor(vi)]), -1)
                    yv = Yt[t.as_tensor(vi)]
                    rec["val_nll"] = float(-lo[t.arange(len(vi)), yv].mean())
                    rec["val_recovery"] = float((lo.argmax(-1) == yv).float().mean())
            self.history.append(rec)
            if verbose:
                print("    " + "  ".join(f"{k} {v:.4f}" if isinstance(v, float) else
                                          f"{k} {v}" for k, v in rec.items()), flush=True)
        return self

    def logprob(self, X: np.ndarray, batch: int = 200_000) -> np.ndarray:
        t = self.torch
        self.net.eval()
        out = []
        with t.no_grad():
            for s in range(0, len(X), batch):
                z = (X[s:s + batch] - self.mu) / self.sd
                out.append(t.log_softmax(self.net(t.as_tensor(z)), -1).numpy())
        return np.concatenate(out) if out else np.zeros((0, NAA), np.float32)

    def save(self, path):
        self.torch.save({"sd": self.net.state_dict(), "mu": self.mu, "s": self.sd,
                         "bg": self.bg, "d_in": self.d_in, "width": self.width,
                         "depth": self.depth, "hist": getattr(self, "history", [])},
                        path)

    def load(self, path):
        z = self.torch.load(path, weights_only=False)
        self.net.load_state_dict(z["sd"])
        self.mu, self.sd, self.bg = z["mu"], z["s"], z["bg"]
        self.history = z.get("hist", [])
        return self


def model_path(variant: str) -> str:
    return os.path.join(MODEL_DIR, f"{variant}.pt")


def load_model(variant: str) -> InvFoldNet:
    """Reconstruct the architecture the checkpoint was saved with, then load it."""
    import torch
    torch.set_num_threads(int(os.environ.get("NT", "2")))
    z = torch.load(model_path(variant), weights_only=False)
    m = InvFoldNet(int(z["d_in"]), width=int(z.get("width", 512)),
                   depth=int(z.get("depth", 3)))
    return m.load(model_path(variant))


#: variant -> (uses neighbour identities, training residues).
#: `geom200k` is a DATA-SCALING probe, not a candidate system. Finding 2 measured that
#: 30x more `prots/` data made the distance predictor MONOTONICALLY worse, because a
#: fragment inside a folded protein has a more compact conformational distribution than an
#: isolated peptide. The central claim of this sprint is that inverse folding should
#: SURVIVE that shift, since local geometry -> residue propensity is close to universal.
#: That claim is testable rather than assertable: train the same model on 200 K residues
#: and on 1.5 M and see which way selected RMSD moves. Scoring it costs almost nothing
#: because it rides along in the same selection pass.
VARIANTS = {
    "geom": (False, 1_500_000),
    "nbid": (True, 1_500_000),
    "geom200k": (False, 200_000),
}


def stage_train(variants=None, epochs: int = 12, width: int = 512):
    os.makedirs(MODEL_DIR, exist_ok=True)
    names = (list(VARIANTS) if variants is None
             else (variants.split(",") if isinstance(variants, str) else list(variants)))
    out = {}
    for v in names:
        if os.path.exists(model_path(v)):
            print(f"{v}: exists")
            continue
        nbid, max_res = VARIANTS[v]
        t0 = time.time()
        X, Y = build_xy(nbid=nbid, max_res=max_res)
        tb = time.time() - t0
        print(f"{v}: X {X.shape} in {tb:.0f}s", flush=True)
        t0 = time.time()
        m = InvFoldNet(X.shape[1], width=int(width)).fit(X, Y, epochs=int(epochs))
        m.save(model_path(v))
        out[v] = {"n_res": int(len(X)), "d": int(X.shape[1]), "nbid": nbid,
                  "build_s": round(tb, 1), "train_s": round(time.time() - t0, 1),
                  "history": m.history}
        del X, Y
        pth = os.path.join(HERE, "invfold_train.json")
        prev = json.load(open(pth)) if os.path.exists(pth) else {}
        prev.update(out)
        json.dump(prev, open(pth, "w"), indent=1)


# ===================================================================== scoring
class SelfConsistency:
    """Score candidate structures by how well they explain a fixed target sequence.

    DEPLOYABLE. Reads the candidate geometry and the target sequence, nothing else.
    """

    def __init__(self, sequence: str, model: InvFoldNet, nbid: bool = False):
        self.seq = sequence
        self.model = model
        self.nbid = nbid
        self.labels = np.array([AA_INDEX.get(a, 0) for a in sequence])
        self.known = np.array([a in AA_INDEX for a in sequence])

    def logprobs(self, coords, phi=None, psi=None) -> np.ndarray:
        X = features(coords, phi, psi,
                     nb_ids=self.labels if self.nbid else None)
        B, n, D = X.shape
        return self.model.logprob(X.reshape(-1, D)).reshape(B, n, NAA)

    def score(self, coords, phi=None, psi=None, pmi: bool = False) -> np.ndarray:
        """``(B,)`` negative mean log-likelihood of the target sequence. LOWER IS BETTER.

        `pmi=True` subtracts the training-set background ``log P(s_i)``, which removes the
        residue-frequency term. That term is IDENTICAL for every candidate of a given
        target, so it cannot change a within-target ranking; it is reported because it
        changes the scale on which targets are compared and is the natural normalisation.
        """
        lp = self.logprobs(coords, phi, psi)
        n = lp.shape[1]
        take = lp[:, np.arange(n), self.labels]
        if pmi:
            take = take - self.model.bg[self.labels][None, :]
        m = self.known
        return -(take[:, m].mean(1) if m.any() else take.mean(1))

    def score_centred(self, coords, phi=None, psi=None) -> np.ndarray:
        """Per-residue ENTROPY-CENTRED variant: subtract each position's own mean log-prob.

        The raw score rewards a candidate whose predicted distributions happen to be
        confident, independently of whether they are confident about the RIGHT residue.
        Centring on ``mean_a log P(a | geometry_i)`` removes that. Unlike the background
        (PMI) subtraction, this term DOES vary between candidates, so it is a genuinely
        different ranking and not a constant shift.
        """
        lp = self.logprobs(coords, phi, psi)
        n = lp.shape[1]
        take = lp[:, np.arange(n), self.labels] - lp.mean(-1)
        m = self.known
        return -(take[:, m].mean(1) if m.any() else take.mean(1))

    def recovery(self, coords, phi=None, psi=None) -> np.ndarray:
        lp = self.logprobs(coords, phi, psi)
        return (lp.argmax(-1) == self.labels[None]).mean(1)


_STAGES = {}


# ===================================================================== evaluation
# Everything below is the EVALUATION HARNESS. It deliberately reads natives, for two
# purposes only: the oracle native-percentile diagnostic (is the objective correctly
# specified at all?) and reporting the CA-RMSD of a selection that was made without them.
# No deployable path imports any of it.
import distogram as dgm                                              # noqa: E402
from s5.lib import B62, encode, kabsch_rmsd_batch, windows_full      # noqa: E402
from s7 import debias as s7d                                         # noqa: E402

K_POOL = int(os.environ.get("K", "500"))
BAND = 1.5
SEL_OUT = os.path.join(HERE, os.environ.get("SEL_OUT", "invfold_select.json"))

_FRAG_MEMO = {}


def _fold_frags(fold):
    if fold not in _FRAG_MEMO:
        _FRAG_MEMO.clear()          # one fold's fragments at a time: ~6k entries each
        _FRAG_MEMO[fold] = list(dgm._fold_fragments(fold, 5))
    return _FRAG_MEMO[fold]


def build_pool(target, folds, K=None):
    """The 126-target instrument's pool for one target, WITH torsions.

    Reproduces `s7.audit.build_target` exactly -- same library (out-of-fold peptides +
    `_fold_fragments`), same window enumeration order, same BLOSUM key, same stable
    argsort -- but through `windows_full`, which additionally carries each window's
    (phi, psi). `ca_raw` is returned so the reproduction can be checked against the
    cached `rr` in `s7/debias_cache/`.
    """
    K = K_POOL if K is None else int(K)
    fold = folds[target.seq]
    peps = [q for q in pdb.load() if folds[q.seq] != fold and q.seq != target.seq]
    CA, S, PHI, PSI, _ = windows_full(peps + _fold_frags(fold), target.n)
    sim = B62[S, encode(target.seq)[None, :]].sum(1)
    idx = np.argsort(-sim, kind="stable")[:K]
    phi, psi = PHI[idx].astype(float), PSI[idx].astype(float)
    # PROVENANCE, second instance, and it bites the decisive experiment specifically.
    # `geo.extract_torsions` fills phi[0] and psi[n-1] with DEFAULT_PHI/DEFAULT_PSI
    # because they are genuinely undefined without the previous C and the next N. A
    # peptide-DB native therefore carries PLACEHOLDERS at both ends, while a window sliced
    # out of a longer parent chain carries the parent's REAL torsions there. Those two
    # values are direct inputs to the torsion feature block, so scoring the native against
    # the pool would partly compare "placeholder vs real end torsion" rather than
    # conformation. Forcing the candidates to the same convention removes it, and costs
    # nothing geometrically: `build_backbone_batch` never uses phi[0] for chain
    # propagation and uses psi[n-1] only to place the terminal O, which is not a feature.
    phi[:, 0] = geo.DEFAULT_PHI
    psi[:, -1] = geo.DEFAULT_PSI
    return {"fold": int(fold), "idx": idx, "sim": sim[idx],
            "ca_raw": CA[idx], "phi": phi, "psi": psi, "n_windows": int(len(CA))}


def _rk(x):
    o = np.argsort(np.argsort(np.asarray(x, float))).astype(float)
    return (o - o.mean()) / max(o.std(), 1e-12)


def partial_rg(score, rms, rg):
    """Spearman(score, rmsd) with radius of gyration partialled out.

    MANDATORY here. Finding 10 recorded that Amber's apparent in-band skill was carried by
    compactness, and finding 1 that the generator-agreement lead went +0.594 -> -0.07 under
    a partial correlation. A scorer that merely likes compact structures gets RMSD
    correlation for free. Identical formula to `s7/dist_baseline.partial_rg`.
    """
    a, r, g = _rk(score), _rk(rms), _rk(rg)
    p_ar, p_ag, p_rg = float(np.mean(a * r)), float(np.mean(a * g)), float(np.mean(r * g))
    den = np.sqrt(max((1 - p_ag ** 2) * (1 - p_rg ** 2), 1e-12))
    return float((p_ar - p_ag * p_rg) / den)


def _spearman(a, b):
    from scipy.stats import spearmanr
    if len(a) < 4:
        return None
    v = spearmanr(a, b).statistic
    return None if not np.isfinite(v) else float(v)


def _arm_stats(s, s_nat, rr, rg, band):
    d = {"sel": float(rr[int(np.argmin(s))]),
         "rho": _spearman(s, rr),
         "rho_band": _spearman(s[band], rr[band]) if band.sum() > 4 else None,
         "rho_partial_rg": partial_rg(s, rr, rg),
         "rho_rg": _spearman(s, rg)}
    if s_nat is not None:
        d["native_pct"] = float(100.0 * (s < s_nat).mean())
        d["native_is_argmin"] = bool(s_nat < s.min())
    return d


def measure_target(target, folds, scorers, K=None, shortlist=(50, 100)):
    """All arms on ONE target's pool. Every arm sees byte-identical candidates."""
    t0 = time.time()
    P = build_pool(target, folds, K)
    phi, psi = P["phi"], P["psi"]
    coords = geo.build_backbone_batch(phi, psi)
    Cca = coords["CA"]

    # ---- ORACLE DIAGNOSTIC region: natives read deliberately, for reporting only ----
    rr = kabsch_rmsd_batch(Cca, target.ca)                  # rebuilt candidates vs native
    rr_raw = kabsch_rmsd_batch(P["ca_raw"], target.ca)      # cache-reproduction check
    rg = np.sqrt(((Cca - Cca.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    nphi = np.asarray(target.phi, float)
    npsi = np.asarray(target.psi, float)
    nat = geo.build_backbone(nphi, npsi)
    nat_b = {k: np.asarray(v)[None] for k, v in nat.items()}
    # ---- end oracle region: rr/rg/nat never enter a ranking ----

    band = rr <= rr.min() + BAND
    rec = {"pdb": target.pdb, "n": int(target.n), "fold": P["fold"], "seq": target.seq,
           "K": int(len(rr)), "n_windows": P["n_windows"],
           "pool_best": float(rr.min()), "pool_mean": float(rr.mean()),
           "pool_best_raw": float(rr_raw.min()), "pool_mean_raw": float(rr_raw.mean()),
           "band_n": int(band.sum()), "arms": {}}

    S = {}
    for name, fn in scorers.items():
        s = np.asarray(fn(target, coords, phi, psi, P), float)
        s_nat = float(np.asarray(
            fn(target, nat_b, nphi[None], npsi[None], P), float)[0])
        S[name] = s
        rec["arms"][name] = _arm_stats(s, s_nat, rr, rg, band)
        if hasattr(fn, "recovery"):
            # The honest held-out sequence-recovery number: this target is excluded from
            # training by sequence identity, and it is a short isolated peptide, which is
            # exactly the distribution the score is deployed on. Reported for the native
            # and for the pool so the gap between them is visible.
            rv = np.asarray(fn.recovery(target, coords, phi, psi, P), float)
            rec["arms"][name]["pool_recovery"] = float(rv.mean())
            rec["arms"][name]["best_recovery"] = float(rv[int(np.argmin(rr))])
            rec["arms"][name]["native_recovery"] = float(np.asarray(
                fn.recovery(target, nat_b, nphi[None], npsi[None], P), float)[0])

    def _z(x):
        return (x - x.mean()) / max(x.std(), 1e-12)

    if "dist" in S:
        # Combine only the two PRIMARY inverse-folding arms. The centred variants and the
        # shuffled control are diagnostics; crossing every one of them with two weights
        # and two shortlists would put ~30 arms in a table read for a winner, and this
        # project has already promoted a proxy once (finding 3, best-MAE arm ranked 16/28).
        for iv in [k for k in ("iv_geom", "iv_nbid") if k in S]:
            for w in (0.5, 1.0):
                rec["arms"]["comb_%s_w%s" % (iv, w)] = _arm_stats(
                    _z(S["dist"]) + w * _z(S[iv]), None, rr, rg, band)
            for M in shortlist:
                sub = np.argsort(S["dist"], kind="stable")[:M]
                pick = sub[int(np.argmin(S[iv][sub]))]
                rec["arms"]["short%d_%s" % (M, iv)] = {
                    "sel": float(rr[pick]),
                    "shortlist_ceiling": float(rr[sub].min()),
                    "shortlist_dist_sel": float(rr[sub[0]]),
                }
    rec["wall"] = round(time.time() - t0, 1)
    return rec


def make_scorers(variants=("geom", "nbid", "geom200k"), with_dist=True):
    """DEPLOYABLE scorers: (target, coords, phi, psi, pool) -> (B,) lower-is-better.

    `target` is read ONLY for its sequence and its identity fold (which selects the
    correct leave-fold-out distogram). No native coordinate, distance or torsion is read
    by any scorer here.
    """
    sc = {}
    models = {v: load_model(v) for v in variants}
    cache = {}

    def _sc(v, seq):
        key = (v, seq)
        if key not in cache:
            cache.clear()
            cache[key] = SelfConsistency(seq, models[v], nbid=VARIANTS[v][0])
        return cache[key]

    for v in variants:
        def f(target, coords, phi, psi, P, v=v):
            return _sc(v, target.seq).score(coords, phi, psi)

        def rec(target, coords, phi, psi, P, v=v):
            return _sc(v, target.seq).recovery(coords, phi, psi)
        f.recovery = rec
        sc["iv_" + v] = f

        def fc(target, coords, phi, psi, P, v=v):
            return _sc(v, target.seq).score_centred(coords, phi, psi)
        sc["iv_%s_ctr" % v] = fc

    # ---- CONTROL, and the most informative single arm in the table ----
    # Finding 3's `gain 0` and finding 8's typicality reference both established that a
    # scorer can look competent while carrying no target-specific information. Here the
    # analogue is exact: score against a fixed random PERMUTATION of the target's own
    # sequence. Same model, same candidates, same amino-acid composition, same length --
    # only the assignment of residues to positions is destroyed. If this ranks as well as
    # the true sequence, the score measures "is this backbone protein-like", not "does
    # this backbone explain THIS sequence", and the whole hypothesis is dead.
    shuf_cache = {}

    def fs(target, coords, phi, psi, P):
        if target.seq not in shuf_cache:
            shuf_cache.clear()
            rng = np.random.default_rng(abs(hash(target.seq)) % (2 ** 31))
            perm = "".join(np.array(list(target.seq))[rng.permutation(len(target.seq))])
            shuf_cache[target.seq] = SelfConsistency(perm, models[variants[0]],
                                                     nbid=VARIANTS[variants[0]][0])
        return shuf_cache[target.seq].score(coords, phi, psi)
    sc["iv_shuffled_ctrl"] = fs

    if with_dist:
        dcache = {}

        def fd(target, coords, phi, psi, P):
            if target.seq not in dcache:
                dcache.clear()
                dcache[target.seq] = dgm.Distogram.for_target(target.seq, fold=P["fold"])
            return np.asarray(dcache[target.seq].score(coords["CA"]), float)
        sc["dist"] = fd
    return sc


def stage_select(limit=None, variants="geom,nbid,geom200k"):
    """The decisive experiment AND the paired selection table, in one pass.

    Both are computed on the same candidates, so the native-percentile result and the
    selection result cannot disagree because of pool construction -- the exact error
    finding 10 documents.
    """
    vs = tuple(variants.split(","))
    targets = s7d.tuning_targets()
    if limit:
        targets = targets[:int(limit)]
    folds = pdb.folds(5)
    scorers = make_scorers(vs)
    rows = json.load(open(SEL_OUT)) if os.path.exists(SEL_OUT) else []
    done = {r["pdb"] for r in rows}
    todo = sorted([p for p in targets if p.pdb not in done],
                  key=lambda p: (folds[p.seq], p.pdb))     # fold-major: memo loads 5x
    print("%d tuning targets, %d to do, K=%d" % (len(targets), len(todo), K_POOL),
          flush=True)
    for k, p in enumerate(todo):
        try:
            rec = measure_target(p, folds, scorers)
        except Exception as exc:                                     # noqa: BLE001
            rec = {"pdb": p.pdb, "error": "%s: %s" % (type(exc).__name__, exc)}
        rows.append(rec)
        tmp = SEL_OUT + ".tmp"
        json.dump(rows, open(tmp, "w"), indent=1)
        os.replace(tmp, SEL_OUT)
        if "error" in rec:
            print("[%3d] %-6s ERROR %s" % (k + 1, p.pdb, rec["error"][:80]), flush=True)
        else:
            a = rec["arms"]
            bits = " ".join(
                "%s sel %5.2f pct %5.1f" % (v, a[v]["sel"], a[v]["native_pct"])
                for v in scorers)
            print("[%3d] %-6s n=%2d pool %5.2f | %s (%.0fs)"
                  % (k + 1, p.pdb, rec["n"], rec["pool_best"], bits, rec["wall"]),
                  flush=True)


def stage_final(which="dev", variants="geom,nbid,geom200k"):
    """The ONE final evaluation, run at the end and never iterated on.

    Protocol: everything is chosen on the 126-target tuning instrument, which is
    cluster-disjoint from both of these. Touching them more than once, or choosing an arm
    by their result, would invalidate them -- so this writes to its own file and the arm
    to report must be named before it runs.
    """
    global SEL_OUT
    targets = pdb.dev_set(24) if which == "dev" else pdb.benchmark()
    SEL_OUT = os.path.join(HERE, "invfold_final_%s.json" % which)
    folds = pdb.folds(5)
    scorers = make_scorers(tuple(variants.split(",")))
    rows = json.load(open(SEL_OUT)) if os.path.exists(SEL_OUT) else []
    done = {r["pdb"] for r in rows}
    todo = sorted([p for p in targets if p.pdb not in done],
                  key=lambda p: (folds[p.seq], p.pdb))
    print("FINAL EVAL on %s: %d targets, %d to do" % (which, len(targets), len(todo)),
          flush=True)
    for k, p in enumerate(todo):
        try:
            rec = measure_target(p, folds, scorers)
        except Exception as exc:                                     # noqa: BLE001
            rec = {"pdb": p.pdb, "error": "%s: %s" % (type(exc).__name__, exc)}
        rows.append(rec)
        tmp = SEL_OUT + ".tmp"
        json.dump(rows, open(tmp, "w"), indent=1)
        os.replace(tmp, SEL_OUT)
        print("[%3d] %-6s %s" % (k + 1, p.pdb,
                                 rec.get("error", "ok pool %.2f" % rec.get("pool_best", 0))),
              flush=True)


_STAGES["select"] = stage_select
_STAGES["final"] = stage_final




# ===================================================================== report
#: The two published benchmarks this experiment exists to be compared against.
#: `s7/FINDINGS.md` finding 5 (distance objective, 126 targets) and finding 10 FINAL
#: (all-atom Amber, 70 paired targets).
PUBLISHED = {
    "distance objective (s7 finding 5, n=126)":
        {"argmin": "3/126", "mean_pct": 36.8, "median_pct": 32.8},
    "all-atom Amber (s7 finding 10 FINAL, n=70)":
        {"argmin": "0/70", "mean_pct": 54.3, "median_pct": 61.5},
}


def _paired(a, b):
    """Paired difference a - b with a 95% t interval and win/loss counts."""
    from scipy import stats
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = a - b
    n = len(d)
    m = float(d.mean())
    se = float(d.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
    h = float(stats.t.ppf(0.975, n - 1) * se) if n > 1 else float("nan")
    t_p = float(stats.ttest_rel(a, b).pvalue) if n > 1 else float("nan")
    return {"n": n, "mean_diff": m, "ci95": [m - h, m + h], "p": t_p,
            "wins_a": int((d < 0).sum()), "wins_b": int((d > 0).sum()),
            "ties": int((d == 0).sum())}


def _agg(rows, arm, key):
    return [r["arms"][arm][key] for r in rows
            if arm in r["arms"] and r["arms"][arm].get(key) is not None]


def _common(rows, arms, key):
    """Rows where EVERY arm has `key`, so every comparison is on the same targets."""
    return [r for r in rows
            if all(arm in r["arms"] and r["arms"][arm].get(key) is not None
                   for arm in arms)]


def stage_report(path=None):
    rows = [r for r in json.load(open(path or SEL_OUT)) if "arms" in r]
    if not rows:
        print("no rows")
        return
    arms = sorted({a for r in rows for a in r["arms"]})
    base = [a for a in arms if not a.startswith(("comb_", "short"))]
    out = {"n_targets": len(rows), "K": rows[0]["K"],
           "pool_best": float(np.mean([r["pool_best"] for r in rows])),
           "pool_mean": float(np.mean([r["pool_mean"] for r in rows])),
           "published_benchmarks": PUBLISHED}

    # ---------------- the decisive experiment ----------------
    print("=" * 78)
    print("NATIVE-PERCENTILE TEST -- is the objective correctly specified?")
    print("=" * 78)
    print("%-22s %6s  %8s  %8s  %10s" % ("objective", "n", "argmin", "mean pct",
                                         "median pct"))
    nat = {}
    for a in base:
        pct = _agg(rows, a, "native_pct")
        am = _agg(rows, a, "native_is_argmin")
        if not pct:
            continue
        nat[a] = {"n": len(pct), "argmin": int(sum(am)), "mean_pct": float(np.mean(pct)),
                  "median_pct": float(np.median(pct)),
                  "pct_below_10": float(np.mean(np.asarray(pct) < 10) * 100),
                  "pct_below_5": float(np.mean(np.asarray(pct) < 5) * 100)}
        print("%-22s %6d  %4d/%-3d  %8.1f  %10.1f"
              % (a, len(pct), sum(am), len(am), np.mean(pct), np.median(pct)))
    for k, v in PUBLISHED.items():
        print("%-22s %6s  %8s  %8.1f  %10.1f"
              % ("[published] " + k.split(" (")[0][:9], v["argmin"].split("/")[1],
                 v["argmin"], v["mean_pct"], v["median_pct"]))
    out["native"] = nat

    # Held-out sequence recovery. NOT the training val_recovery, which is split by residue
    # over overlapping windows of the same chain and reads far too high. These targets are
    # excluded from training by sequence identity and are short isolated peptides, which is
    # the distribution the score is deployed on.
    recs = {a: {k: (float(np.mean(_agg(rows, a, k))) if _agg(rows, a, k) else None)
                for k in ("native_recovery", "best_recovery", "pool_recovery")}
            for a in base}
    recs = {a: v for a, v in recs.items() if v["native_recovery"] is not None}
    if recs:
        print("")
        print("HELD-OUT SEQUENCE RECOVERY (random = 5.0%, single most frequent aa "
              "~9%):")
        print("%-24s %10s %14s %12s" %
              ("arm", "native", "pool best", "pool mean"))
        for a in sorted(recs, key=lambda x: -recs[x]["native_recovery"]):
            v = recs[a]
            print("%-24s %9.1f%% %13.1f%% %11.1f%%"
                  % (a, 100 * v["native_recovery"], 100 * v["best_recovery"],
                     100 * v["pool_recovery"]))
        out["recovery"] = recs

    # paired native-percentile differences against the distogram on identical pools
    if "dist" in nat:
        out["native_paired_vs_dist"] = {}
        print("")
        print("PAIRED native percentile vs the distogram, identical pools "
              "(negative = ranks the native better):")
        print("%-24s %5s %9s %-22s %8s %9s"
              % ("arm", "n", "d(pct)", "95% CI", "p", "W/L/T"))
        for a in sorted(nat, key=lambda x: nat[x]["mean_pct"]):
            if a == "dist":
                continue
            cm = _common(rows, [a, "dist"], "native_pct")
            pr = _paired([r["arms"][a]["native_pct"] for r in cm],
                         [r["arms"]["dist"]["native_pct"] for r in cm])
            out["native_paired_vs_dist"][a] = pr
            print("%-24s %5d %+9.2f [%+7.2f,%+7.2f] %8.4f  %3d/%3d/%d"
                  % (a, pr["n"], pr["mean_diff"], pr["ci95"][0], pr["ci95"][1], pr["p"],
                     pr["wins_a"], pr["wins_b"], pr["ties"]))

    # ---------------- selection ----------------
    print("\n" + "=" * 78)
    print("SELECTED CA-RMSD (PRIMARY) -- paired, identical candidates")
    print("=" * 78)
    sel = {}
    for a in arms:
        v = _agg(rows, a, "sel")
        if not v:
            continue
        v = np.asarray(v)
        sel[a] = {"n": len(v), "mean": float(v.mean()), "median": float(np.median(v)),
                  "sd": float(v.std(ddof=1)),
                  "se": float(v.std(ddof=1) / np.sqrt(len(v))),
                  "best": float(v.min()), "worst": float(v.max()),
                  "frac_below_2": float((v < 2.0).mean()),
                  "frac_below_1_5": float((v < 1.5).mean())}
    order = sorted(sel, key=lambda a: sel[a]["mean"])
    print("%-24s %5s %7s %7s %7s %7s %7s %7s" %
          ("arm", "n", "mean", "median", "SE", "<2.0A", "<1.5A", "worst"))
    print("%-24s %5d %7.3f %7.3f %7s %7s %7s %7s" %
          ("pool best (ceiling)", len(rows), out["pool_best"],
           float(np.median([r["pool_best"] for r in rows])), "", "", "", ""))
    for a in order:
        s = sel[a]
        print("%-24s %5d %7.3f %7.3f %7.3f %7.1f%% %6.1f%% %7.3f"
              % (a, s["n"], s["mean"], s["median"], s["se"],
                 100 * s["frac_below_2"], 100 * s["frac_below_1_5"], s["worst"]))
    out["selected"] = sel

    if "dist" in sel:
        print("\nPAIRED vs the distogram on identical candidates "
              "(negative favours the arm):")
        print("%-24s %5s %9s %-22s %8s %9s" %
              ("arm", "n", "d(sel)", "95% CI", "p", "W/L/T"))
        out["selected_paired_vs_dist"] = {}
        for a in order:
            if a == "dist":
                continue
            cm = _common(rows, [a, "dist"], "sel")
            if not cm:
                continue
            p = _paired([r["arms"][a]["sel"] for r in cm],
                        [r["arms"]["dist"]["sel"] for r in cm])
            out["selected_paired_vs_dist"][a] = p
            print("%-24s %5d %+9.3f [%+7.3f,%+7.3f] %8.4f  %3d/%3d/%d"
                  % (a, p["n"], p["mean_diff"], p["ci95"][0], p["ci95"][1], p["p"],
                     p["wins_a"], p["wins_b"], p["ties"]))

    # ---------------- correlations, with the compactness control ----------------
    print("\n" + "=" * 78)
    print("RANKING CORRELATIONS -- raw, in-band, and with rg PARTIALLED OUT")
    print("=" * 78)
    print("%-24s %5s %8s %9s %13s %8s" %
          ("arm", "n", "rho", "rho_band", "rho_partial_rg", "rho_rg"))
    corr = {}
    for a in base + [x for x in arms if x.startswith("comb_")]:
        d = {k: (float(np.mean(_agg(rows, a, k))) if _agg(rows, a, k) else None)
             for k in ("rho", "rho_band", "rho_partial_rg", "rho_rg")}
        if d["rho"] is None:
            continue
        corr[a] = d
        print("%-24s %5d %+8.3f %+9.3f %+13.3f %+8.3f"
              % (a, len(_agg(rows, a, "rho")), d["rho"],
                 d["rho_band"] if d["rho_band"] is not None else float("nan"),
                 d["rho_partial_rg"], d["rho_rg"]))
    out["correlations"] = corr
    if "dist" in corr:
        out["corr_paired_vs_dist"] = {}
        print("")
        print("PAIRED vs the distogram (positive favours the arm). The "
              "rg-partialled column")
        print("is the one that matters: it is where finding 10 took Amber apart.")
        print("%-24s %-16s %9s %-22s %9s"
              % ("arm", "statistic", "d", "95% CI", "W/L"))
        for a in sorted(corr):
            if a == "dist":
                continue
            for key in ("rho_band", "rho_partial_rg"):
                cm = _common(rows, [a, "dist"], key)
                if len(cm) <= 2:
                    continue
                pr = _paired([r["arms"][a][key] for r in cm],
                             [r["arms"]["dist"][key] for r in cm])
                out["corr_paired_vs_dist"].setdefault(a, {})[key] = pr
                # NOTE the sign convention flip against the selected-RMSD table above:
                # for RMSD lower is better, for a correlation HIGHER is better. So the
                # difference is printed as-is (arm minus distogram) and a WIN is d > 0.
                print("%-24s %-16s %+9.3f [%+7.3f,%+7.3f] %4d/%-4d"
                      % (a, key, pr["mean_diff"], pr["ci95"][0], pr["ci95"][1],
                         pr["wins_b"], pr["wins_a"]))

    # ---------------- shortlist ceilings ----------------
    sh = [a for a in arms if a.startswith("short")]
    if sh:
        print("\n" + "=" * 78)
        print("SHORTLIST ARMS -- reported WITH the shortlist's own ceiling, so a gain "
              "cannot\nbe confused with a ceiling artefact")
        print("=" * 78)
        print("%-24s %5s %8s %10s %12s" %
              ("arm", "n", "sel", "ceiling", "dist alone"))
        out["shortlist"] = {}
        for a in sorted(sh):
            v = _agg(rows, a, "sel")
            c = _agg(rows, a, "shortlist_ceiling")
            ds = _agg(rows, a, "shortlist_dist_sel")
            out["shortlist"][a] = {"n": len(v), "sel": float(np.mean(v)),
                                   "ceiling": float(np.mean(c)),
                                   "dist_alone": float(np.mean(ds))}
            print("%-24s %5d %8.3f %10.3f %12.3f"
                  % (a, len(v), np.mean(v), np.mean(c), np.mean(ds)))

    # ---------------- named contrasts ----------------
    # The four questions the arm table cannot answer by inspection, each as an explicit
    # paired test. The shuffled-sequence pair is the decisive one.
    CONTRASTS = [
        ("iv_nbid", "iv_shuffled_ctrl", "true sequence vs SHUFFLED (nbid)"),
        ("iv_geom", "iv_shuffled_ctrl", "true sequence vs SHUFFLED (geom)"),
        ("iv_nbid", "iv_geom", "neighbour identity vs geometry only"),
        ("iv_geom", "iv_geom200k", "1.5M vs 200k training residues"),
        ("comb_iv_nbid_w0.5", "iv_nbid", "combination vs inverse folding alone"),
        ("comb_iv_nbid_w0.5", "dist", "combination vs the distogram"),
    ]
    print("")
    print("=" * 78)
    print("NAMED CONTRASTS (negative favours the FIRST arm on RMSD/percentile)")
    print("=" * 78)
    out["contrasts"] = {}
    for key, unit in (("sel", "A"), ("native_pct", "pct")):
        print("-- %s --" % key)
        for a, b, lbl in CONTRASTS:
            cm = _common(rows, [a, b], key)
            if not cm:
                continue
            pr = _paired([r["arms"][a][key] for r in cm],
                         [r["arms"][b][key] for r in cm])
            out["contrasts"].setdefault(key, {})["%s_vs_%s" % (a, b)] = pr
            print("%-38s %+8.3f %-3s [%+7.3f,%+7.3f] p=%.4f  W/L %3d/%3d"
                  % (lbl, pr["mean_diff"], unit, pr["ci95"][0], pr["ci95"][1], pr["p"],
                     pr["wins_a"], pr["wins_b"]))

    # Does the native explain its OWN sequence better than the pool does? This is the
    # hypothesis in its purest form, with no ranking and no RMSD in it at all.
    print("-- held-out recovery, native minus pool mean (POSITIVE favours the "
          "hypothesis) --")
    out["contrasts"]["recovery_native_minus_pool"] = {}
    for a in [x for x in base if _agg(rows, x, "native_recovery")]:
        cm = _common(rows, [a], "native_recovery")
        pr = _paired([r["arms"][a]["native_recovery"] for r in cm],
                     [r["arms"][a]["pool_recovery"] for r in cm])
        out["contrasts"]["recovery_native_minus_pool"][a] = pr
        print("%-38s %+8.2f pts [%+7.2f,%+7.2f] p=%.4f  W/L %3d/%3d"
              % (a, 100 * pr["mean_diff"], 100 * pr["ci95"][0], 100 * pr["ci95"][1],
                 pr["p"], pr["wins_b"], pr["wins_a"]))

    p = os.path.join(HERE, "invfold_report.json")
    json.dump(out, open(p, "w"), indent=1)
    print("\nwrote %s  (n=%d targets, SE on selected = %.3f A)"
          % (p, len(rows), sel.get("dist", {}).get("se", float("nan"))))
    return out


def stage_chirality(n=200):
    """The chirality-sensitivity check, as a REPORTABLE number rather than a test.

    Reflects real pool candidates (not random torsions) and reports how far the
    inverse-folding score moves against a distance score that provably cannot move at all.
    ORACLE-FREE: no native is read.
    """
    folds = pdb.folds(5)
    n = int(n)
    targets = s7d.tuning_targets()[:8]
    models = {v: load_model(v) for v in VARIANTS if os.path.exists(model_path(v))}
    rows = []
    for t in targets:
        P = build_pool(t, folds, K=n)
        c = geo.build_backbone_batch(P["phi"], P["psi"])
        cm = {k: v * np.array([1.0, 1.0, -1.0]) for k, v in c.items()}
        D = np.linalg.norm(c["CA"][:, :, None] - c["CA"][:, None], axis=-1)
        Dm = np.linalg.norm(cm["CA"][:, :, None] - cm["CA"][:, None], axis=-1)
        d = dgm.Distogram.for_target(t.seq, fold=P["fold"])
        sd = np.asarray(d.score(c["CA"]), float)
        sdm = np.asarray(d.score(cm["CA"]), float)
        rec = {"pdb": t.pdb, "n": int(t.n), "K": int(len(sd)),
               "dist_matrix_max_change": float(np.abs(D - Dm).max()),
               "distogram_max_change": float(np.abs(sd - sdm).max())}
        for v, m in models.items():
            sc = SelfConsistency(t.seq, m, nbid=VARIANTS[v][0])
            a = sc.score(c, P["phi"], P["psi"])
            b = sc.score(cm, P["phi"], P["psi"])
            rec["iv_%s_mean_abs_change" % v] = float(np.abs(a - b).mean())
            rec["iv_%s_score_sd" % v] = float(a.std())
            rec["iv_%s_rank_agreement" % v] = float(
                np.mean(np.argsort(np.argsort(a)) == np.argsort(np.argsort(b))))
            rec["iv_%s_mirror_worse_frac" % v] = float((b > a).mean())
        rows.append(rec)
        print(rec, flush=True)
    agg = {k: float(np.mean([r[k] for r in rows]))
           for k in rows[0] if isinstance(rows[0][k], float)}
    out = {"per_target": rows, "mean": agg}
    json.dump(out, open(os.path.join(HERE, "invfold_chirality.json"), "w"), indent=1)
    print("\n=== MEAN ===")
    for k, v in agg.items():
        print("  %-32s %.6g" % (k, v))
    return out


def stage_verify(n_targets=3, K=200):
    """Prove the pool reproduction is the SAME instrument as `s7/debias_cache/`.

    `build_pool` re-derives the pool through `windows_full` instead of `windows_of` so it
    can carry torsions. If the window ORDER differed by even one element the whole
    comparison would be against a different instrument, so this checks the raw-window
    RMSDs against the cached ones element by element.
    """
    n_targets, K = int(n_targets), int(K)
    folds = pdb.folds(5)
    ok, bad = [], []
    for t in s7d.tuning_targets():
        p = os.path.join(HERE, "..", "s7", "debias_cache", "%s.npz" % t.pdb)
        if not os.path.exists(p):
            continue
        z = np.load(p, allow_pickle=True)
        P = build_pool(t, folds, K=K)
        mine = kabsch_rmsd_batch(P["ca_raw"], t.ca)
        theirs = np.asarray(z["rr"], float)[:K]
        e = float(np.abs(mine - theirs).max())
        (ok if e < 1e-4 else bad).append({"pdb": t.pdb, "max_abs_rr_diff": e,
                                          "n_windows": P["n_windows"],
                                          "cached_n_windows": int(z["n_windows"])})
        print("%-6s max|rr diff| %.3e  windows %d vs cached %d"
              % (t.pdb, e, P["n_windows"], int(z["n_windows"])), flush=True)
        if len(ok) + len(bad) >= n_targets:
            break
    out = {"matched": ok, "mismatched": bad, "K": K}
    json.dump(out, open(os.path.join(HERE, "invfold_verify.json"), "w"), indent=1)
    print("\n%d matched, %d mismatched" % (len(ok), len(bad)))
    return out


_STAGES["report"] = stage_report
_STAGES["chirality"] = stage_chirality
_STAGES["verify"] = stage_verify


# ===================================================================== main
def main(argv):
    cmd = argv[1] if len(argv) > 1 else "report"
    if cmd == "leak":
        stage_leak()
    elif cmd == "scan":
        stage_scan()
    elif cmd == "windows":
        stage_windows()
    elif cmd == "train":
        stage_train()
    elif cmd in _STAGES:
        _STAGES[cmd](*argv[2:])
    else:
        print(f"unknown stage {cmd}; have leak scan windows train "
              + " ".join(sorted(_STAGES)))


if __name__ == "__main__":
    main(sys.argv)
