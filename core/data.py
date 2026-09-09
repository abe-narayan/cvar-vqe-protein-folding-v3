"""Sequences, structure databases, encodings and bitstring representations.

Consolidates `peptide_db.py` (434), `fragment_db.py` (145), `representations.py` (633) and
`esm_features.py` (148).

The sharpest hazard in this repository is the sequence alphabet, so it is dealt with first.
Two different 20-letter orders were in use and they agree on exactly THREE letters (A, S,
T):

    s7.audit.encode, s5.lib, s6.library, s8.generate   "ARNDCQEGHILKMFPSTWYV"
    s9.evo.encode, s8.invfold, invfold.py              "ACDEFGHIKLMNPQRSTVWY"

The cached window banks in ``s8/generate_univ/*.npz`` are written in the FIRST order.
`s10/forensics.stage_identity` decoded them with the SECOND, and therefore measured, for
every one of 126 targets, the identity between the target and a systematically scrambled
sequence -- 1CEK's own exact-sequence window decoded as ``ALSYMMAGAYQMM`` instead of
``AISVLLAQAVFLL`` and scored 0.308 instead of 1.000. The whole leakage audit was void.

Two structural defences, not one:

1. `ALPHABET` is the canonical order, pinned by `tests/test_data.py` against the deposited
   banks themselves -- the test decodes cached windows and requires them to be substrings
   of their own targets, which the wrong alphabet fails.
2. `decode` takes the alphabet as a REQUIRED argument and `CodedBank` carries it in the
   file next to the codes. There is no way to decode without saying which alphabet, and no
   way to save codes without saving the alphabet with them. The 2026 failure needed a
   default; there is now no default to be wrong about.

The identity convention. `identity` is Needleman-Wunsch match count normalised by the
LONGER sequence, applied member-to-member. This is leaky at the member level for long
library entries: 573 library members across 58 targets contain a target at >= 0.6 while
passing the 0.6 filter, because 1CEK's 13 residues sitting verbatim inside the 25-residue
1A11 scores 13/25 = 0.520. Results on record depend on this convention, so it is preserved
EXACTLY and pinned by a test. `containment` is the same alignment normalised by the SHORTER
sequence -- the quantity a window-level filter should actually use -- and is provided under
its own name so the two can never be confused for each other again.

Retrieval tie-breaking is pinned. `s7/audit.py` sorts candidates with a STABLE argsort and
`s5/inband.py` does not. Measured over the 126 cached universes, the stable sort reproduces
the deposited ``order`` on every target exactly, while the unstable sort moves as many as
47 of 500 pool members. `top_k` is stable, always.
"""
from __future__ import annotations

import glob
import json
import math
import os
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import cache
from . import geometry as geo

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 1. the alphabet
#: THE canonical residue order. Everything in `core` encodes and decodes in this order,
#: and it is the order the cached window banks on disk are already written in.
ALPHABET = "ARNDCQEGHILKMFPSTWYV"
#: The other order that exists in the repository. Kept ONLY so a test can demonstrate that
#: it is the wrong one; never used to encode or decode anything.
ALPHABET_ALT = "ACDEFGHIKLMNPQRSTVWY"
_IDX = {a: i for i, a in enumerate(ALPHABET)}
N_AA = len(ALPHABET)


def encode(seq: str) -> np.ndarray:
    """``(n,)`` int8 codes in `ALPHABET` order. Unknown letters map to 0, as they always
    have -- changing that would move pool membership."""
    return np.array([_IDX.get(c, 0) for c in seq], dtype=np.int64)


def decode(codes: np.ndarray, alphabet: str) -> str:
    """Codes -> letters. ``alphabet`` is REQUIRED and has no default, on purpose: the one
    time it had a default, an audit decoded 126 targets' banks with the wrong one."""
    if len(alphabet) != N_AA:
        raise ValueError(f"alphabet must have {N_AA} letters, got {len(alphabet)}")
    return "".join(alphabet[int(c)] for c in np.asarray(codes).ravel())


class CodedBank:
    """Integer-encoded sequences that carry their own alphabet.

    Saving codes without the alphabet is what made the void audit possible, so this object
    refuses to be constructed or written without one, and `load` refuses a file that lacks
    it.
    """

    def __init__(self, codes: np.ndarray, alphabet: str = ALPHABET):
        if len(alphabet) != N_AA:
            raise ValueError("alphabet must have 20 letters")
        self.codes = np.asarray(codes)
        self.alphabet = alphabet

    @classmethod
    def of(cls, sequences: Sequence[str]) -> "CodedBank":
        return cls(np.stack([encode(s) for s in sequences]))

    def decode(self, row: int = 0) -> str:
        return decode(self.codes[row], self.alphabet)

    def decode_all(self) -> List[str]:
        return [decode(r, self.alphabet) for r in np.atleast_2d(self.codes)]

    def save(self, path: str) -> None:
        np.savez_compressed(path, codes=self.codes,
                            alphabet=np.array(self.alphabet))

    @classmethod
    def load(cls, path: str) -> "CodedBank":
        z = np.load(path, allow_pickle=False)
        if "alphabet" not in z.files:
            raise ValueError(
                f"{path} stores codes without an alphabet. Refusing to guess: the two "
                f"orders in this repository agree on three letters out of twenty.")
        return cls(z["codes"], str(z["alphabet"]))

#: BLOSUM62 in `ALPHABET` order, upper triangle flattened.
_B62_ROWS = """4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
4 2 -3 1 0 -3 -2 -1 -3 -1 3
4 -2 2 0 -3 -2 -1 -2 -1 1
5 -1 -3 -1 0 -1 -3 -2 -2
5 0 -2 -1 -1 -1 -1 1
6 -4 -2 -2 1 3 -1
7 -1 -1 -4 -3 -2
4 1 -3 -2 -2
5 -2 -2 0
11 2 -3
7 -1
4"""


def _blosum62() -> np.ndarray:
    M = np.zeros((N_AA, N_AA))
    for i, row in enumerate(r.split() for r in _B62_ROWS.strip().split("\n")):
        for k, v in enumerate(row):
            M[i, i + k] = M[i + k, i] = float(v)
    return M

BLOSUM62 = _blosum62()


def blosum_similarity(query_codes: np.ndarray, bank_codes: np.ndarray) -> np.ndarray:
    """``(W,)`` summed BLOSUM62 similarity of every row of ``bank_codes`` to the query.

    One gather over a (20, 20) table rather than a Python loop over windows: the bank is
    ``(W, n)`` and the query ``(n,)``, so this is ``BLOSUM62[query, bank].sum(-1)`` with the
    query broadcast.
    """
    q = np.asarray(query_codes).ravel()
    B = np.atleast_2d(np.asarray(bank_codes))
    return BLOSUM62[q[None, :], B].sum(-1)


def top_k(scores: np.ndarray, k: int, largest: bool = True) -> np.ndarray:
    """The pinned retrieval ordering: a STABLE argsort, so ties resolve by bank position.

    `s7/audit.py` used ``kind="stable"`` and `s5/inband.py` did not. Over the 126 cached
    window universes the stable sort reproduces the deposited ``order`` array on every
    target with zero mismatches; the unstable sort changes up to 47 of 500 pool members,
    which is a different pool and therefore a different result. There is no defensible
    reason to leave that to the sorting algorithm, so it is fixed here.
    """
    s = np.asarray(scores, float)
    order = np.argsort(-s if largest else s, kind="stable")
    return order[:k] if k is not None else order


# 2. Alignment and identity
def identity(a: str, b: str, gap: float = -1.0) -> float:
    """Needleman-Wunsch identity normalised by the LONGER sequence.

    The repository's convention, preserved exactly -- results on record depend on it. It
    is leaky at the member level (see the module docstring); `containment` is the
    normalisation a window-level filter should use, and it lives under its own name so the
    two cannot be mistaken for one another.
    """
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return 0.0
    return _nw_matches(a, [b])[0] / max(n, m)


def containment(a: str, b: str, gap: float = -1.0) -> float:
    """The same alignment normalised by the SHORTER sequence.

    For a retrieval pool the comparable unit is the window, so the filter should be applied
    window-to-target (equal lengths, where the two normalisations coincide) or equivalently
    member-to-target normalised by the shorter sequence. Normalising by the longer one lets
    a 26-residue library peptide carry a 9-residue target verbatim at 9/26 = 0.346.
    """
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return 0.0
    return _nw_matches(a, [b])[0] / min(n, m)

#: Below this many sequences the per-cell numpy overhead of the batched DP exceeds the
#: pure-Python inner loop it replaces (the batched form does ``n*m`` numpy calls whatever
#: the batch size). Measured crossover on this database is around 24; the two paths are
#: asserted bit-identical in `tests/test_data.py`, so this is purely a speed switch.
_BATCH_MIN = 24


def _nw_matches_scalar(a: str, b: str, gap: float = -1.0) -> float:
    """The original per-pair Needleman-Wunsch, kept verbatim as the reference form."""
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return 0.0
    prev = [(gap * j, 0) for j in range(m + 1)]
    for i in range(1, n + 1):
        cur = [(gap * i, 0)]
        ai = a[i - 1]
        for j in range(1, m + 1):
            match = 1.0 if ai == b[j - 1] else 0.0
            best = (prev[j - 1][0] + match, prev[j - 1][1] + int(match > 0))
            up = (prev[j][0] + gap, prev[j][1])
            left = (cur[j - 1][0] + gap, cur[j - 1][1])
            if up[0] > best[0]:
                best = up
            if left[0] > best[0]:
                best = left
            cur.append(best)
        prev = cur
    return float(prev[m][1])


def _nw_matches(a: str, bs: Sequence[str], gap: float = -1.0) -> np.ndarray:
    """Match counts of ``a`` against every sequence in ``bs``, batched over ``bs``.

    The DP is sequential in both indices, but completely independent ACROSS pairs, so the
    inner cell update runs on the whole batch at once: the ``n*m`` Python cell updates
    become ``n*m`` numpy operations on ``B``-vectors. 83 us/pair scalar -> 3.4 us/pair at
    B = 787, and the whole 309k-pair clustering pass drops from 25.7 s to 3.4 s.

    Ragged lengths are handled by RIGHT-padding with a byte that matches no residue and
    reading each row out at its own column ``len(bs[k])``. That is exact, not an
    approximation: cell ``(i, j)`` depends only on cells with column ``<= j``, so columns
    past a sequence's own length cannot influence its answer.

    The tie-break order is the scalar routine's, exactly: ``diag`` wins, then ``up`` on a
    strict improvement, then ``left`` on a strict improvement. Changing it would change
    which alignment is reported and therefore the match count.
    """
    n = len(a)
    B = len(bs)
    out = np.zeros(B)
    if n == 0 or B == 0:
        return out
    if B < _BATCH_MIN:
        for k, b in enumerate(bs):
            out[k] = _nw_matches_scalar(a, b, gap)
        return out
    lens = np.array([len(b) for b in bs])
    m = int(lens.max())
    if m == 0:
        return out
    ai = np.frombuffer(a.encode("latin-1"), dtype=np.uint8)
    bj = np.zeros((B, m), dtype=np.uint8)          # 0 is not a residue letter
    for k, b in enumerate(bs):
        if b:
            bj[k, :len(b)] = np.frombuffer(b.encode("latin-1"), dtype=np.uint8)

    prev_s = np.repeat((gap * np.arange(m + 1, dtype=float))[None, :], B, 0)
    prev_c = np.zeros((B, m + 1))
    cur_s = np.empty_like(prev_s)
    cur_c = np.empty_like(prev_c)
    for i in range(1, n + 1):
        cur_s[:, 0] = gap * i
        cur_c[:, 0] = 0.0
        match = (bj == ai[i - 1]).astype(float)                         # (B, m)
        ps, pc = prev_s, prev_c
        for j in range(1, m + 1):
            mv = match[:, j - 1]
            best_s = ps[:, j - 1] + mv
            best_c = pc[:, j - 1] + mv
            up_s = ps[:, j] + gap
            take = up_s > best_s
            best_s = np.where(take, up_s, best_s)
            best_c = np.where(take, pc[:, j], best_c)
            le_s = cur_s[:, j - 1] + gap
            take = le_s > best_s
            cur_s[:, j] = np.where(take, le_s, best_s)
            cur_c[:, j] = np.where(take, cur_c[:, j - 1], best_c)
        prev_s, cur_s = cur_s, prev_s
        prev_c, cur_c = cur_c, prev_c
    return prev_c[np.arange(B), lens]


def identity_many(a: str, bs: Sequence[str], gap: float = -1.0) -> np.ndarray:
    """`identity` of ``a`` against many sequences at once. Same numbers, batched."""
    if not len(bs):
        return np.zeros(0)
    n = len(a)
    lens = np.array([len(b) for b in bs], float)
    denom = np.maximum(lens, n)
    return np.where(denom > 0, _nw_matches(a, bs, gap) / np.maximum(denom, 1), 0.0)


def _composition(s: str) -> Dict[str, int]:
    c: Dict[str, int] = {}
    for ch in s:
        c[ch] = c.get(ch, 0) + 1
    return c


def max_possible_identity(a: str, b: str, ca: Optional[Dict[str, int]] = None,
                          cb: Optional[Dict[str, int]] = None) -> float:
    """An ADMISSIBLE upper bound on `identity(a, b)`.

    No alignment can match a residue type more often than it occurs in the scarcer of the
    two sequences, so the multiset intersection bounds the match count from above. This is
    a sound prefilter -- it can only over-estimate -- unlike the 3-mer Jaccard test it
    replaced, whose premise ("identity above 0.6 forces a shared 3-mer") is false because
    matches may alternate with mismatches. That unsound filter let one benchmark target
    train on a 0.70-identity homolog. It skips ~98% of pairs.
    """
    ca = ca or _composition(a)
    cb = cb or _composition(b)
    shared = sum(min(v, cb.get(ch, 0)) for ch, v in ca.items())
    return shared / max(len(a), len(b))


def composition_matrix(seqs: Sequence[str]) -> Tuple[np.ndarray, np.ndarray]:
    """``(counts (N, 20) int16, lengths (N,) int32)`` -- the prefilter in array form."""
    C = np.zeros((len(seqs), N_AA), dtype=np.int16)
    for k, s in enumerate(seqs):
        if s:
            np.add.at(C[k], encode(s), 1)
    return C, np.array([len(s) for s in seqs], dtype=np.int32)


def max_possible_identity_many(a: str, C: np.ndarray, lens: np.ndarray) -> np.ndarray:
    """`max_possible_identity` of ``a`` against a whole `composition_matrix`, vectorised.

    The bound is the multiset intersection over the alphabet, which is
    ``minimum(count_a, count_b).sum(1)`` -- one numpy reduction instead of 20 dictionary
    lookups per pair. Over the 309k clustering pairs that is 0.68 s of Python replaced by
    9 ms of numpy, and the value is identical because both compute the same integer.
    """
    ca = np.zeros(N_AA, dtype=np.int16)
    if a:
        np.add.at(ca, encode(a), 1)
    shared = np.minimum(ca[None, :], C).sum(1)
    return shared / np.maximum(np.maximum(lens, len(a)), 1)


def _kmers(s: str, k: int = 3) -> set:
    return {s[i:i + k] for i in range(len(s) - k + 1)}

# 3. The peptide database
DIRS = (os.path.join(BASE, "pdbs"), os.path.join(BASE, "pdbs_ext"))
PEPTIDE_CACHE = os.path.join(BASE, "peptide_db.npz")

MIN_LEN, MAX_LEN = 8, 26
#: Backbone RMSD between a native and the structure rebuilt from its own torsions by the
#: ideal-geometry builder. Above this the deposited geometry carries something the
#: representation cannot express (cis bonds, distorted valence), so the entry is not
#: evidence about anything the search can reach.
REBUILD_TOL = 1.5
IDENTITY_THRESHOLD = 0.6


@dataclass
class Peptide:
    pdb: str
    seq: str
    ca: np.ndarray            # (n, 3)
    phi: np.ndarray           # (n,) radians
    psi: np.ndarray
    rebuild: float

    @property
    def n(self) -> int:
        return len(self.seq)


def _scan() -> List[Dict]:
    seen: Dict[str, Dict] = {}
    for d in DIRS:
        for path in sorted(glob.glob(os.path.join(d, "*.pdb"))):
            pdbid = os.path.basename(path)[:-4].upper()
            try:
                seq, coords, phi, psi = geo.native_coords_from_pdb(path)
            except Exception:
                continue
            n = len(seq)
            if not (MIN_LEN <= n <= MAX_LEN):
                continue
            if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(psi))):
                continue
            ca = np.asarray(coords["CA"], float)
            if len(ca) != n:
                continue
            step = np.linalg.norm(np.diff(ca, axis=0), axis=1)
            if step.min() < 3.5 or step.max() > 4.1:
                continue          # chain break or a non-standard linkage
            try:
                built = geo.build_backbone(phi, psi)
                keys = [k for k in ("N", "CA", "C") if k in coords and k in built]
                mob = np.vstack([built[k] for k in keys])
                tgt = np.vstack([coords[k] for k in keys])
                rb = float(geo.rmsd(geo.kabsch_superpose(mob, tgt), tgt))
            except Exception:
                continue
            if not np.isfinite(rb) or rb > REBUILD_TOL:
                continue
            prev = seen.get(seq)
            if prev is None or rb < prev["rebuild"]:
                seen[seq] = dict(pdb=pdbid, seq=seq, ca=ca, phi=phi, psi=psi, rebuild=rb)
    return sorted(seen.values(), key=lambda r: r["pdb"])


def build_peptides(force: bool = False) -> List[Peptide]:
    if not force and os.path.exists(PEPTIDE_CACHE):
        z = np.load(PEPTIDE_CACHE, allow_pickle=True)
        return [Peptide(**r) for r in z["records"]]
    recs = _scan()
    out = [Peptide(**r) for r in recs]
    np.savez_compressed(PEPTIDE_CACHE, records=np.array(recs, dtype=object))
    return out


@lru_cache(maxsize=1)
def load() -> Tuple[Peptide, ...]:
    return tuple(build_peptides())


@lru_cache(maxsize=1)
def _seq_index() -> Dict[str, int]:
    return {p.seq: k for k, p in enumerate(load())}


def by_pdb(pdbid: str) -> Optional[Peptide]:
    return _by_pdb().get(pdbid.upper())


@lru_cache(maxsize=1)
def _by_pdb() -> Dict[str, Peptide]:
    return {p.pdb: p for p in load()}


def holdout(target_seq: str, threshold: float = IDENTITY_THRESHOLD) -> Tuple[Peptide, ...]:
    """Database with every entry too similar to `target_seq` removed.

    EXACT: every surviving database entry is aligned against the target. The 3-mer
    prefilter this used to carry was unsound, and one benchmark target was provably trained
    on a 0.70-identity homolog because of it.
    """
    target = target_seq.strip().upper()
    entries = load()
    C, lens = _db_composition()
    bound = max_possible_identity_many(target, C, lens)
    seqs = [p.seq for p in entries]
    out, maybe, maybe_i = [], [], []
    for k, s in enumerate(seqs):
        if s == target:
            continue
        if bound[k] < threshold:
            out.append(k)
        else:
            maybe.append(s)
            maybe_i.append(k)
    if maybe:
        for k, v in zip(maybe_i, identity_many(target, maybe)):
            if v < threshold:
                out.append(k)
    return tuple(entries[k] for k in sorted(out))


@lru_cache(maxsize=1)
def _db_composition() -> Tuple[np.ndarray, np.ndarray]:
    return composition_matrix([p.seq for p in load()])

# ------------------------------------------------------------------ clusters and folds
_CLUSTER_CACHE = os.path.join(BASE, "peptide_clusters.json")
_FOLD_CACHE = os.path.join(BASE, "peptide_folds.json")


def clusters(threshold: float = IDENTITY_THRESHOLD) -> Dict[str, int]:
    """Single-linkage identity clusters over the whole database, sequence -> cluster id.

    EXACT: every pair that survives the admissible composition bound is aligned. The
    previous k-mer prefilter was unsound and 1C9A, a reported benchmark target, was trained
    on 2N0G at 0.70 identity because of it.
    """
    if os.path.exists(_CLUSTER_CACHE):
        with open(_CLUSTER_CACHE) as f:
            saved = json.load(f)
        if saved.get("threshold") == threshold:
            return saved["assign"]
    seqs = [p.seq for p in load()]
    C, lens = composition_matrix(seqs)
    parent = list(range(len(seqs)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a in range(len(seqs)):
        bound = max_possible_identity_many(seqs[a], C, lens)
        bound[:a + 1] = 0.0
        cand_i = [int(b) for b in np.flatnonzero(bound >= threshold) if find(a) != find(b)]
        if not cand_i:
            continue
        for b, v in zip(cand_i, identity_many(seqs[a], [seqs[b] for b in cand_i])):
            if v >= threshold and find(a) != find(b):
                parent[find(b)] = find(a)
    roots: Dict[int, int] = {}
    assign: Dict[str, int] = {}
    for i, s in enumerate(seqs):
        assign[s] = roots.setdefault(find(i), len(roots))
    with open(_CLUSTER_CACHE, "w") as f:
        json.dump({"threshold": threshold, "assign": assign}, f)
    return assign


def folds(n_folds: int = 5, seed: int = 0) -> Dict[str, int]:
    """Sequence -> fold index, assigning whole identity clusters to folds. PINNED.

    The assignment is written to `peptide_folds.json` on first use and read back
    thereafter, and that pinning is load-bearing rather than a cache. Fold indices are a
    shuffle of cluster ids, so anything that renumbers clusters -- adding a database entry,
    or replacing an unsound prefilter with a sound one -- silently reassigns every sequence,
    the trained fold models on disk would then be excluding the wrong fold, and a target
    could be scored by a model that had trained on it with no error raised anywhere. Delete
    this file only together with every model in `distogram_models/`, `pairnet_models/` and
    `torsion_models/`.
    """
    if os.path.exists(_FOLD_CACHE):
        with open(_FOLD_CACHE) as f:
            saved = json.load(f)
        if saved.get("n_folds") == n_folds and saved.get("seed") == seed:
            assign = saved["assign"]
            if not [p.seq for p in load() if p.seq not in assign]:
                return assign
    assign = clusters()
    cl = sorted(set(assign.values()))
    rng = np.random.default_rng(seed)
    rng.shuffle(cl)
    fold_of = {c: i % n_folds for i, c in enumerate(cl)}
    out = {s: fold_of[c] for s, c in assign.items()}
    with open(_FOLD_CACHE, "w") as f:
        json.dump({"n_folds": n_folds, "seed": seed, "assign": out}, f)
    return out

# ------------------------------------------------------------------ target sets
MANIFEST = os.path.join(BASE, "results", "benchmark_manifest.json")
MONOMER_MANIFEST = os.path.join(BASE, "results", "monomer_manifest.json")
BENCHMARK_N = 60
BENCHMARK_MIN_LEN = 9
BENCHMARK_MAX_LEN = 16


def benchmark_set(min_len: int = 9, max_len: int = 16,
                  pdbs: Optional[Sequence[str]] = None) -> List[Peptide]:
    """Targets used for reporting. Defaults to the repo's own `pdbs/` entries in range."""
    if pdbs:
        keep_seq = {p.seq for p in load() if p.pdb in set(pdbs)}
    else:
        keep_seq = set()
        for path in sorted(glob.glob(os.path.join(BASE, "pdbs", "*.pdb"))):
            try:
                seq, _, _, _ = geo.parse_pdb(path)
            except Exception:
                continue
            keep_seq.add(seq)
    # Deduplication may have replaced a local entry with a better-rebuilding deposit of the
    # SAME sequence (1J4M -> 1K43, 5AWL -> 2RVD). Selecting by sequence keeps the target.
    return [p for p in load() if p.seq in keep_seq and min_len <= p.n <= max_len]


def dev_set(n: int = 24, min_len: int = 9, max_len: int = 16, seed: int = 7,
            _seq_only: bool = False) -> List[Peptide]:
    """Hyperparameter development targets, CLUSTER-disjoint from `benchmark`.

    An earlier version excluded only exact benchmark sequences, which left four development
    peptides sharing an identity cluster with six reported targets. Hyperparameters fitted
    on that set measurably transferred to its homologs and not to the rest: a controlled
    re-run helped the six contaminated targets by 0.135 A and hurt eight controls by
    0.289 A, a difference-in-differences of 0.424 A.

    ``_seq_only=True`` is the internal acyclic form: `benchmark` needs a dev set to exclude
    while it is being defined. Benchmark membership is fixed by the internal form and does
    not depend on the stronger filter, so it is stable.
    """
    cl = clusters()
    if _seq_only:
        bench_seq = {p.seq for p in benchmark_set(min_len, max_len)}
        pool = [p for p in load()
                if p.seq not in bench_seq and min_len <= p.n <= max_len]
    else:
        bench_cl = {cl[p.seq] for p in benchmark()}
        pool = [p for p in load()
                if cl[p.seq] not in bench_cl and min_len <= p.n <= max_len]
    by_cluster: Dict[int, Peptide] = {}
    for p in sorted(pool, key=lambda x: x.pdb):
        by_cluster.setdefault(cl[p.seq], p)
    reps = sorted(by_cluster.values(), key=lambda x: x.pdb)
    idx = np.random.default_rng(seed).permutation(len(reps))[:n]
    return [reps[i] for i in sorted(idx)]


def benchmark(n: int = BENCHMARK_N, min_len: int = BENCHMARK_MIN_LEN,
              max_len: int = BENCHMARK_MAX_LEN, seed: int = 11) -> List[Peptide]:
    """The reported benchmark, FROZEN to `results/benchmark_manifest.json` once written.

    Membership is generated once and then pinned. This is not a convenience: the selection
    depends on identity-cluster ids, and correcting the clustering silently swapped 13 of
    the 60 targets. A benchmark whose membership moves when an internal detail changes
    cannot support a before/after claim, and re-running the baseline is not the fix.
    """
    if os.path.exists(MANIFEST):
        with open(MANIFEST) as f:
            want = [t["pdb"] for t in json.load(f)["targets"]]
        by_id = _by_pdb()
        missing = [w for w in want if w not in by_id]
        if missing:
            raise KeyError(f"manifest names targets absent from the database: {missing}")
        return [by_id[w] for w in want]
    core = benchmark_set(min_len, max_len)
    dev = {p.seq for p in dev_set(24, min_len=9, max_len=16, _seq_only=True)}
    cl = clusters()
    taken_seq = {p.seq for p in core}
    taken_cl = {cl[p.seq] for p in core}
    pool = [p for p in load()
            if min_len <= p.n <= max_len and p.seq not in taken_seq
            and p.seq not in dev and cl[p.seq] not in taken_cl]
    by_cluster: Dict[int, Peptide] = {}
    for p in sorted(pool, key=lambda x: x.pdb):
        by_cluster.setdefault(cl[p.seq], p)
    reps = sorted(by_cluster.values(), key=lambda x: x.pdb)
    order = np.random.default_rng(seed).permutation(len(reps))
    out = list(core)
    for i in order:
        if len(out) >= n:
            break
        out.append(reps[i])
    return sorted(out, key=lambda p: (p.n, p.pdb))


def monomer_benchmark(aqueous_only: bool = False) -> List[Peptide]:
    """The 52-target monomer benchmark, or the 36-target strictly-aqueous subset.

    The 60-target set is dominated by chains that are not folded on their own -- fibril
    segments (7.22 A), micelle- and membrane-bound peptides (5.27 A), partner-bound
    (4.03 A) -- where the reference conformation is held by inter-chain packing or a
    detergent, and a single-chain intramolecular energy has no access to it. Selection is
    on the EXPERIMENT (deposited header text, model count, ensemble spread, builder
    reproducibility, cluster uniqueness), never on any method's performance.
    """
    with open(MONOMER_MANIFEST) as f:
        man = json.load(f)
    by_id = _by_pdb()
    out = []
    for t in man["targets"]:
        if aqueous_only and not t.get("aqueous", True):
            continue
        if t["pdb"] in by_id:
            out.append(by_id[t["pdb"]])
    return out

# 4. FRAGMENTS
PROT_DIR = os.path.join(BASE, "prots")
FRAGMENT_CACHE = os.path.join(BASE, "fragment_db.npz")
FRAGMENT_CACHE_LARGE = os.path.join(BASE, "fragment_db_large.npz")
FRAG_LENGTHS = tuple(range(9, 21))
FRAG_STRIDE = 5
FRAG_REBUILD_TOL = 1.0
#: Set FRAG_LARGE=1 to make every consumer read the larger extraction. A switch rather than
#: a silent default, so a result can always be attributed to a fragment set.
USE_LARGE = os.environ.get("FRAG_LARGE", "") == "1"


def _extract_fragments(path: str) -> List[Dict]:
    try:
        seq, coords, phi, psi = geo.native_coords_from_pdb(path)
    except Exception:
        return []
    ca = np.asarray(coords["CA"], float)
    n = len(seq)
    if n < min(FRAG_LENGTHS) or len(ca) != n:
        return []
    step = np.linalg.norm(np.diff(ca, axis=0), axis=1)
    ok = (step > 3.5) & (step < 4.1)
    stem = os.path.basename(path)[:-4].upper()
    out = []
    for L in FRAG_LENGTHS:
        for s in range(1, n - L, FRAG_STRIDE):    # skip residue 0: phi is undefined
            e = s + L
            if e >= n:
                break
            if not ok[s:e - 1].all():
                continue
            ph, ps = phi[s:e], psi[s:e]
            if not (np.all(np.isfinite(ph)) and np.all(np.isfinite(ps))):
                continue
            sub = seq[s:e]
            if "X" in sub:
                continue
            out.append(dict(pdb=f"{stem}_{s}", seq=sub, ca=ca[s:e].copy(),
                            phi=ph.copy(), psi=ps.copy(), rebuild=0.0))
    return out


def build_fragments(max_per_protein: int = 6, n_target: int = 6000, seed: int = 0,
                    force: bool = False, path: Optional[str] = None) -> List[Peptide]:
    path = path or FRAGMENT_CACHE
    if not force and os.path.exists(path):
        z = np.load(path, allow_pickle=True)
        return [Peptide(**r) for r in z["records"]]
    rng = np.random.default_rng(seed)
    seen: Dict[str, Dict] = {}
    paths = sorted(glob.glob(os.path.join(PROT_DIR, "*.pdb")))
    rng.shuffle(paths)
    for p in paths:
        frags = _extract_fragments(p)
        if not frags:
            continue
        rng.shuffle(frags)
        kept = 0
        for f in frags:
            if kept >= max_per_protein:
                break
            if f["seq"] in seen:
                continue
            # Only fragments the ideal-geometry builder can express are evidence about
            # anything the search can build.
            try:
                built = geo.build_backbone(f["phi"], f["psi"])
                if geo.rmsd(geo.kabsch_superpose(built["CA"], f["ca"]),
                            f["ca"]) > FRAG_REBUILD_TOL:
                    continue
            except Exception:
                continue
            seen[f["seq"]] = f
            kept += 1
        if len(seen) >= n_target:
            break
    recs = list(seen.values())
    np.savez_compressed(path, records=np.array(recs, dtype=object))
    return [Peptide(**r) for r in recs]


@lru_cache(maxsize=2)
def load_fragments(large: Optional[bool] = None) -> Tuple[Peptide, ...]:
    if large is None:
        large = USE_LARGE
    return tuple(build_fragments(path=FRAGMENT_CACHE_LARGE if large else FRAGMENT_CACHE))

_FOLD_FRAG_VERSION = 1


def fold_fragments(fold: int, n_folds: int = 5,
                   threshold: float = IDENTITY_THRESHOLD,
                   large: Optional[bool] = None) -> List[Peptide]:
    """Fragments safe to train on for `fold`: none above the identity threshold to any
    peptide the fold holds out.

    Cached. It was 3.0 s per fold, recomputed on every `train_fold` call for both model
    families -- 30 s of pure re-derivation per sweep -- and it is a pure function of
    (fold, n_folds, threshold, fragment set).
    """
    if large is None:
        large = USE_LARGE
    frags = load_fragments(large)
    held = [p.seq for p in load() if folds(n_folds)[p.seq] == fold]
    params = dict(fold=int(fold), n_folds=int(n_folds), threshold=float(threshold),
                  large=bool(large), n_frag=len(frags), n_held=len(held),
                  held_digest=_digest(held), frag_digest=_digest([f.seq for f in frags]))

    def compute():
        hk = [_kmers(s) for s in held]
        keep = []
        for k, f in enumerate(frags):
            fk = _kmers(f.seq)
            cand = [s for s, kk in zip(held, hk) if not (fk and kk and not (fk & kk))]
            if cand and (identity_many(f.seq, cand) >= threshold).any():
                continue
            keep.append(k)
        return {"keep": np.array(keep, dtype=np.int64)}

    got = cache.cached("fold_frag", _FOLD_FRAG_VERSION, compute, **params)
    return [frags[k] for k in got["keep"]]


def _digest(items: Sequence[str]) -> str:
    import hashlib
    h = hashlib.blake2b(digest_size=16)
    for s in items:
        h.update(s.encode())
        h.update(b"\x00")
    return h.hexdigest()


# 5. WINDOWS -- the retrieval unit
def windows(pool: Sequence[Peptide], n: int, with_torsions: bool = False):
    """Every contiguous length-``n`` window of every pool member.

    Returns ``(CA (W, n, 3), codes (W, n), src [(member, start)])`` and, with
    ``with_torsions``, ``phi`` and ``psi`` as well. A window carries its parent's torsions
    directly, so a retrieved candidate has a real conformation with nothing to infer:
    ``phi[s]`` of a window starting at ``s`` is defined by the residue before it in the
    parent chain, which exists.

    Codes are in `ALPHABET` order. Wrap them in a `CodedBank` before persisting.
    """
    cas, seqs, phis, psis, src = [], [], [], [], []
    for k, q in enumerate(pool):
        m = len(q.seq)
        if m < n:
            continue
        e = encode(q.seq)
        for s in range(m - n + 1):
            cas.append(q.ca[s:s + n])
            seqs.append(e[s:s + n])
            src.append((k, s))
            if with_torsions:
                phis.append(q.phi[s:s + n])
                psis.append(q.psi[s:s + n])
    if not cas:
        raise ValueError(f"no pool member is at least {n} residues long")
    out = [np.stack(cas), np.stack(seqs), src]
    if with_torsions:
        out += [np.stack(phis).astype(np.float32), np.stack(psis).astype(np.float32)]
    return tuple(out)


def retrieve(query_seq: str, pool: Sequence[Peptide], k: int):
    """Top-``k`` BLOSUM62-similar windows for ``query_seq``. Stable ties -- see `top_k`."""
    n = len(query_seq)
    W, S, src = windows(pool, n)
    sim = blosum_similarity(encode(query_seq), S)
    idx = top_k(sim, k)
    return W[idx], S[idx], sim[idx], [src[i] for i in idx]

# 6. ESM-2 FEATURES
#: The full bank. ~1.5 GB, an object-array npz: `np.load` materialises ALL of it, and doing
#: that in a process that then keeps running has repeatedly taken this box to 94-96%. It is
#: never opened from a long-lived process here -- `warm` extracts the needed subset in a
#: throwaway subprocess (`cache.extract_subset`).
ESM_CACHE = os.path.join(BASE, "esm_cache.npz")
#: The hot cache: only the sequences inference actually needs, ~21 MB.
ESM_SMALL = os.path.join(BASE, "esm_small.npz")
ESM_PCA = os.path.join(BASE, "esm_pca.npz")
ESM_MODEL = "esm2_t33_650M_UR50D"
N_PCA = 32

_small: Optional[Dict[str, Tuple[np.ndarray, np.ndarray]]] = None
_pca_cache: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None


def _read_npz_pairs(path: str) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    if not os.path.exists(path):
        return {}
    z = np.load(path, allow_pickle=True)
    return {str(k): (v[0], v[1]) for k, v in zip(z["seqs"], z["vals"])}


def _hot() -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    global _small
    if _small is None:
        _small = _read_npz_pairs(ESM_SMALL)
    return _small


def esm_available() -> bool:
    return os.path.exists(ESM_CACHE) or os.path.exists(ESM_SMALL)


def warm(sequences: Sequence[str]) -> int:
    """Make ``sequences`` available from the hot cache, extracting from the 1.5 GB bank in
    a throwaway subprocess. Returns how many were added.

    This is the established pattern and the reason it exists: caching per-target tables
    instead of the whole bank cut one downstream footprint from 2.0 GB to 0.26 GB.
    """
    global _small
    have = _hot()
    todo = [s for s in dict.fromkeys(sequences) if s not in have]
    if not todo:
        return 0
    if not os.path.exists(ESM_CACHE):
        raise FileNotFoundError(ESM_CACHE)
    tmp = os.path.join(BASE, "core_cache", "_esm_subset.npz")
    cache.extract_subset(ESM_CACHE, todo, tmp)
    got = _read_npz_pairs(tmp)
    have.update(got)
    keys = list(have.keys())
    np.savez_compressed(ESM_SMALL, seqs=np.array(keys, dtype=object),
                        vals=np.array([[have[s][0], have[s][1]] for s in keys],
                                      dtype=object))
    os.remove(tmp)
    _small = have
    return len(got)


def esm_compute(sequences: Sequence[str], batch: int = 8, verbose: bool = True) -> None:
    """Run ESM-2 over `sequences` and extend the full bank. Requires torch + fair-esm."""
    import torch
    import esm

    bank = _read_npz_pairs(ESM_CACHE)
    todo = [s for s in dict.fromkeys(sequences) if s not in bank]
    if not todo:
        return
    model, alphabet = getattr(esm.pretrained, ESM_MODEL)()
    model.eval()
    bc = alphabet.get_batch_converter()
    layer = model.num_layers
    with torch.no_grad():
        for start in range(0, len(todo), batch):
            chunk = todo[start:start + batch]
            _, _, toks = bc([(f"p{i}", s) for i, s in enumerate(chunk)])
            out = model(toks, repr_layers=[layer], return_contacts=True)
            reps = out["representations"][layer].numpy()
            cons = out["contacts"].numpy()
            for m, s in enumerate(chunk):
                n = len(s)
                bank[s] = (reps[m, 1:n + 1].astype(np.float32),
                           cons[m, :n, :n].astype(np.float32))
            if verbose:
                print(f"  esm {start + len(chunk)}/{len(todo)}", flush=True)
    seqs = list(bank.keys())
    np.savez_compressed(ESM_CACHE, seqs=np.array(seqs, dtype=object),
                        vals=np.array([[bank[s][0], bank[s][1]] for s in seqs],
                                      dtype=object))


def esm_raw(sequence: str) -> Tuple[np.ndarray, np.ndarray]:
    """``(per-residue representation (n, D), contact probabilities (n, n))``.

    Reads the hot cache; falls back to extracting from the big bank in a subprocess, and
    only then to running the model. Never loads the big bank into this process.
    """
    hot = _hot()
    if sequence in hot:
        return hot[sequence]
    if os.path.exists(ESM_CACHE):
        warm([sequence])
        hot = _hot()
        if sequence in hot:
            return hot[sequence]
    esm_compute([sequence], verbose=False)
    warm([sequence])
    return _hot()[sequence]


def fit_esm_pca(sequences: Sequence[str], n_components: int = N_PCA) -> None:
    X = np.vstack([esm_raw(s)[0] for s in sequences])
    mu = X.mean(0)
    U, S, Vt = np.linalg.svd(X - mu, full_matrices=False)
    np.savez(ESM_PCA, mu=mu, W=Vt[:n_components].T,
             scale=(S[:n_components] / np.sqrt(len(X))))


def _pca() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    global _pca_cache
    if _pca_cache is None:
        z = np.load(ESM_PCA)
        _pca_cache = (z["mu"], z["W"], z["scale"])
    return _pca_cache

_ESM_VERSION = 1


def esm_embed(sequence: str) -> np.ndarray:
    """``(n, N_PCA)`` whitened per-residue embedding.

    A 1280-dimensional embedding against ~13,000 database residues would let a downstream
    tree model memorise identity rather than learn a mapping; the PCA is a fixed,
    sequence-independent projection onto the directions the peptide corpus actually varies
    along. Cached per sequence, because the projection is deterministic and the whole point
    of the hot cache is to never touch the 1.5 GB bank twice.
    """
    def compute():
        mu, W, scale = _pca()
        return {"e": ((esm_raw(sequence)[0] - mu) @ W) / np.maximum(scale, 1e-6)}
    return cache.cached("esm_embed", _ESM_VERSION, compute,
                        seq=sequence, n_pca=N_PCA,
                        pca_mtime=repr(os.path.getmtime(ESM_PCA))
                        if os.path.exists(ESM_PCA) else "none")["e"]


def esm_contacts(sequence: str) -> np.ndarray:
    """ESM-2's own contact head: a logistic regression over symmetrised, APC-corrected
    attention maps trained on real contacts -- a pairwise quantity, which is the shape the
    distance prior needs, and the only part of ESM-2 ever supervised on structure."""
    return esm_raw(sequence)[1]

# 7. REPRESENTATIONS -- bitstring to structure
import sidechains as _sc                                                   # noqa: E402

#: One-letter codes whose chi1 is encoded: the aromatics whose ring this repo can build.
#: Chi1 rotates the ring about CA-CB, the degree of freedom that decides whether a stacked
#: aromatic pair is geometrically reachable -- one qubit per aromatic residue, against ten
#: to move a 10-mer from 4 to 8 backbone states for 0.04 A of ceiling.
CHI1_ENCODED = tuple(sorted(geo.THREE_TO_ONE[r] for r in _sc.AROMATIC_RING_RESIDUES))

# Convention relied on by callers that build reference structures: index 0 is the helical
# state of its class and index 1 the extended/beta state, in EVERY library. So
# `bitstring_from_states([0] * n)` is "as helical as this sequence can be" regardless of
# composition. A single global library was applied to every residue before this: glycine
# was offered nothing isoleucine was not, and proline was allowed phi = +60.
GENERAL_4 = [(-63.0, -42.0), (-120.0, 130.0), (-75.0, 150.0), (-90.0, 0.0)]
GENERAL_8 = [(-63.0, -42.0), (-120.0, 130.0), (-75.0, 150.0), (-140.0, 160.0),
             (-60.0, -30.0), (-90.0, 0.0), (-100.0, 120.0), (60.0, 45.0)]
#: Glycine has no C-beta, so the left-handed and mirror regions are genuinely populated --
#: for Gly they are where a large fraction of residues sit, not exotic.
GLY_4 = [(-63.0, -42.0), (-120.0, 130.0), (60.0, 45.0), (90.0, -175.0)]
GLY_8 = [(-63.0, -42.0), (-120.0, 130.0), (60.0, 45.0), (90.0, -175.0),
         (-90.0, 0.0), (75.0, -10.0), (-80.0, 165.0), (105.0, 170.0)]
#: Proline: phi is locked near -63 by the pyrrolidine ring, which cannot open to +60. The
#: two real clusters are alpha-R-like and PPII, resolved more finely instead of offering
#: geometrically impossible states.
PRO_4 = [(-63.0, -35.0), (-63.0, 150.0), (-70.0, 120.0), (-63.0, -10.0)]
PRO_8 = [(-63.0, -35.0), (-63.0, 150.0), (-70.0, 120.0), (-63.0, -10.0),
         (-55.0, 145.0), (-75.0, 160.0), (-60.0, -50.0), (-70.0, 175.0)]
#: The residue N-terminal to a proline: the following ring depletes alpha-R, enhances
#: beta/PPII and opens the otherwise-unused zeta region. Modelling this is why the library
#: has to be sequence-CONTEXT aware, not merely residue-type aware.
PRE_PRO_4 = [(-70.0, -35.0), (-130.0, 135.0), (-75.0, 155.0), (-140.0, 80.0)]
PRE_PRO_8 = [(-70.0, -35.0), (-130.0, 135.0), (-75.0, 155.0), (-140.0, 80.0),
             (-100.0, 125.0), (-155.0, 160.0), (-85.0, 105.0), (-60.0, 140.0)]
#: The exact pre-2026-07-28 flat libraries, preserved verbatim. `run_8state_seed0.py` and
#: `amber_obc2.py` are pinned to a persisted set of bitstring -> energy pairs that
#: `partest/golden_check.py` verifies to max absolute error 0.0; changing the library
#: changes the coordinates a bitstring decodes to and therefore every energy in that set.
#: Do not "clean this up" -- that silently invalidates a 4-hour reference run.
LEGACY_4 = [(-63.0, -42.0), (-120.0, 130.0), (-75.0, 150.0), (60.0, 45.0)]
LEGACY_8 = [(-63.0, -42.0), (-120.0, 130.0), (-75.0, 150.0), (-140.0, 160.0),
            (-60.0, -30.0), (-90.0, 0.0), (60.0, 45.0), (90.0, -10.0)]
LEGACY_LIBRARIES = {4: LEGACY_4, 8: LEGACY_8}

CLASS_GENERAL, CLASS_GLY, CLASS_PRO = "GENERAL", "GLY", "PRO"
CLASS_PRE_PRO, CLASS_LEGACY = "PRE_PRO", "LEGACY"

STATE_LIBRARIES: Dict[int, Dict[str, List[Tuple[float, float]]]] = {
    4: {CLASS_GENERAL: GENERAL_4, CLASS_GLY: GLY_4,
        CLASS_PRO: PRO_4, CLASS_PRE_PRO: PRE_PRO_4},
    8: {CLASS_GENERAL: GENERAL_8, CLASS_GLY: GLY_8,
        CLASS_PRO: PRO_8, CLASS_PRE_PRO: PRE_PRO_8},
}
STATES_4, STATES_8 = GENERAL_4, GENERAL_8


def residue_classes(sequence: Optional[str], n_residues: int) -> List[str]:
    """Per-residue library class. Falls back to GENERAL everywhere without a sequence."""
    if not sequence:
        return [CLASS_GENERAL] * n_residues
    seq = sequence.strip().upper()
    if len(seq) != n_residues:
        raise ValueError(f"sequence length {len(seq)} != n_residues {n_residues}")
    out = []
    for i, aa in enumerate(seq):
        nxt = seq[i + 1] if i + 1 < len(seq) else ""
        out.append(CLASS_GLY if aa == "G" else
                   CLASS_PRO if aa == "P" else
                   CLASS_PRE_PRO if nxt == "P" else CLASS_GENERAL)
    return out


def _bits_needed(k: int) -> int:
    b = 1
    while (1 << b) < k:
        b += 1
    return b


def _ang_diff_deg(a, b):
    return ((a - b + 180.0) % 360.0) - 180.0


class TorsionStateRepresentation:
    """Per-residue (phi, psi) drawn from a discrete, residue-class-specific library."""

    name = "torsion"
    is_lattice = False

    def __init__(self, n_residues: int, n_states: int = 4,
                 sequence: Optional[str] = None, legacy_library: bool = False,
                 chi_bits: bool = True):
        if n_states not in STATE_LIBRARIES:
            raise ValueError(f"n_states must be one of {sorted(STATE_LIBRARIES)}")
        self.n_residues = int(n_residues)
        self.n_states = int(n_states)
        self.sequence = sequence.strip().upper() if sequence else None
        self.legacy_library = bool(legacy_library)
        if self.legacy_library:
            # `legacy_library` means the whole ENCODING is the old one, so chi1 bits are off
            # too: adding them would change n_bits and the bitstring -> structure map.
            self.classes = [CLASS_LEGACY] * self.n_residues
            self.libraries = [LEGACY_LIBRARIES[n_states]] * self.n_residues
            self.chi_bits = False
        else:
            self.classes = residue_classes(self.sequence, self.n_residues)
            libs = STATE_LIBRARIES[n_states]
            self.libraries = [libs[c] for c in self.classes]
            self.chi_bits = bool(chi_bits)

        self.bits_per_residue = _bits_needed(self.n_states)
        self.n_backbone_bits = self.bits_per_residue * self.n_residues
        # Chi bits are appended AFTER the backbone bits so the backbone slot layout is
        # untouched: search moves and the ceiling search address slots by
        # `i * bits_per_residue` and must keep working unchanged.
        self.chi_residues: Tuple[int, ...] = tuple(
            i for i, aa in enumerate(self.sequence or "") if aa in CHI1_ENCODED
        ) if self.chi_bits else tuple()
        self.n_chi_bits = len(self.chi_residues)
        self.n_bits = self.n_backbone_bits + self.n_chi_bits
        self.n_qubits = self.n_bits
        self.offsets = [i * self.bits_per_residue for i in range(self.n_residues)]
        self._rows = np.arange(self.n_residues)
        self._phi = np.radians(np.array([[p for p, _ in lib] for lib in self.libraries]))
        self._psi = np.radians(np.array([[q for _, q in lib] for lib in self.libraries]))
        self._lib_deg = np.array([[[p, q] for p, q in lib] for lib in self.libraries])
        self._ring_spec = tuple(
            (i, geo.ONE_TO_THREE[self.sequence[i]],
             _sc.ring_atom_names(geo.ONE_TO_THREE[self.sequence[i]]))
            for i in self.chi_residues)
        self._rotamers = tuple(_sc.CHI1_ROTAMERS[geo.ONE_TO_THREE[self.sequence[i]]]
                               for i in self.chi_residues)

    # -- unchecked internals; a single validated entry point does all the work once ----
    def _states(self, bitstring: str) -> List[int]:
        b = self.bits_per_residue
        return [int(bitstring[o:o + b], 2) % self.n_states for o in self.offsets]

    def _angles(self, idx) -> Tuple[np.ndarray, np.ndarray]:
        return self._phi[self._rows, idx].copy(), self._psi[self._rows, idx].copy()

    def _chi(self, bitstring: str) -> Dict[int, float]:
        off = self.n_backbone_bits
        return {i: float(rot[int(bitstring[off + k]) % len(rot)])
                for k, ((i, _, _), rot) in enumerate(zip(self._ring_spec, self._rotamers))}

    def _rings(self, coords, chi) -> Dict[int, Dict[str, np.ndarray]]:
        out: Dict[int, Dict[str, np.ndarray]] = {}
        for i, key, names in self._ring_spec:
            if not names:
                continue
            atoms = _sc.build_sidechain(key, coords["N"][i], coords["CA"][i],
                                        coords["C"][i], coords["CB"][i], chi1=chi.get(i))
            out[i] = {nm: atoms[nm] for nm in names if nm in atoms}
        return out

    def _check(self, bitstring: str) -> None:
        if len(bitstring) != self.n_bits:
            raise ValueError(f"bitstring length {len(bitstring)} != n_bits {self.n_bits}")

    # -- decode ---------------------------------------------------------------------
    def state_indices(self, bitstring: str) -> List[int]:
        self._check(bitstring)
        return self._states(bitstring)

    def decode(self, bitstring: str) -> Tuple[np.ndarray, np.ndarray]:
        """bitstring -> (phi, psi) in radians.

        Does NOT overwrite phi[0] and psi[n-1] with constants: `build_backbone` never reads
        phi[0] for any CA position, but `energy_terms.rama_penalty` WAS being evaluated on
        it, so residue 0's torsion penalty was computed from a torsion the structure did not
        have.
        """
        self._check(bitstring)
        return self._angles(self._states(bitstring))

    def build_coords(self, bitstring: str) -> Dict[str, np.ndarray]:
        return geo.build_backbone(*self.decode(bitstring))

    def build_all(self, bitstring: str):
        """``(phi, psi, coords, rings)`` in one pass. The Hamiltonian's entry point:
        validates once, decodes once, builds the backbone once, and builds the rings from
        that same backbone. `rings` is None when chi1 is not encoded."""
        self._check(bitstring)
        phi, psi = self._angles(self._states(bitstring))
        coords = geo.build_backbone(phi, psi)
        rings = self._rings(coords, self._chi(bitstring)) if self._ring_spec else None
        return phi, psi, coords, (rings or None)

    def decode_batch(self, states: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """``(B, n)`` state indices -> two ``(B, n)`` angle arrays. One gather, no loop."""
        s = np.asarray(states, int)
        return self._phi[self._rows[None], s], self._psi[self._rows[None], s]

    def build_batch(self, states: np.ndarray) -> Dict[str, np.ndarray]:
        return geo.build_backbone_batch(*self.decode_batch(states))

    # -- chi1 -----------------------------------------------------------------------
    def chi_indices(self, bitstring: str) -> List[int]:
        self._check(bitstring)
        off = self.n_backbone_bits
        return [int(bitstring[off + k]) for k in range(self.n_chi_bits)]

    def chi1_degrees(self, bitstring: str) -> Dict[int, float]:
        if not self.n_chi_bits:
            return {}
        self._check(bitstring)
        return self._chi(bitstring)

    def build_rings(self, bitstring: str, coords=None) -> Dict[int, Dict[str, np.ndarray]]:
        """Aromatic ring atoms per residue, honouring the encoded chi1. Only the residues
        in `chi_residues` are built -- typically two out of ten. Prefer `build_all` on the
        hot path; it avoids rebuilding the backbone."""
        if not self.n_chi_bits:
            return {}
        self._check(bitstring)
        if coords is None:
            coords = geo.build_backbone(*self._angles(self._states(bitstring)))
        return self._rings(coords, self._chi(bitstring))

    # -- encode ---------------------------------------------------------------------
    def _nearest_state(self, i: int, phi_rad, psi_rad) -> int:
        """Closest library state for residue `i`. ``phi_rad``/``psi_rad`` may be None at a
        terminus, where the torsion is genuinely undefined -- scoring on a placeholder there
        would bias the snapped native."""
        lib = self._lib_deg[i]
        d = np.zeros(len(lib))
        if phi_rad is not None:
            d += _ang_diff_deg(math.degrees(phi_rad), lib[:, 0]) ** 2
        if psi_rad is not None:
            d += _ang_diff_deg(math.degrees(psi_rad), lib[:, 1]) ** 2
        return int(d.argmin())

    def bitstring_from_states(self, states: Sequence[int],
                              chi_states: Optional[Sequence[int]] = None) -> str:
        w = self.bits_per_residue
        if len(states) != self.n_residues:
            raise ValueError(f"got {len(states)} states for {self.n_residues} residues")
        bits = "".join(format(int(s) % self.n_states, f"0{w}b") for s in states)
        if not self.n_chi_bits:
            return bits
        if chi_states is None:
            chi = [0] * self.n_chi_bits
        else:
            chi = list(chi_states)
            # Silently truncating would emit a short bitstring that every downstream
            # `_check` then rejects with a confusing message about n_bits.
            if len(chi) != self.n_chi_bits:
                raise ValueError(f"got {len(chi)} chi states for {self.n_chi_bits} chi "
                                 f"bits (residues {self.chi_residues})")
        return bits + "".join(str(int(c) % 2) for c in chi)

    def angles_to_bits(self, phi, psi, chi_states=None) -> str:
        n = self.n_residues
        return self.bitstring_from_states(
            [self._nearest_state(i, None if i == 0 else phi[i],
                                 None if i == n - 1 else psi[i]) for i in range(n)],
            chi_states=chi_states)

    def native_bitstring(self, native_phi, native_psi, **_) -> str:
        """Native backbone, default chi rotamers. Chi1 is deliberately not snapped from the
        native: only N, CA and C are read from a PDB here, and chi1 does not affect CA
        positions, so the representation ceiling is unaffected either way."""
        return self.angles_to_bits(native_phi, native_psi)

    def random_bitstring(self, rng: np.random.Generator) -> str:
        return self.bitstring_from_states(
            rng.integers(0, self.n_states, size=self.n_residues),
            chi_states=(rng.integers(0, 2, size=self.n_chi_bits)
                        if self.n_chi_bits else None))

    def describe(self) -> Dict[str, object]:
        return {"name": self.name, "n_states": self.n_states,
                "bits_per_residue": self.bits_per_residue, "n_qubits": self.n_qubits,
                "n_backbone_bits": self.n_backbone_bits, "n_chi_bits": self.n_chi_bits,
                "chi_residues": list(self.chi_residues),
                "config_space": (float(self.n_states) ** self.n_residues
                                 * float(2 ** self.n_chi_bits)),
                "residue_classes": list(self.classes),
                "sequence_aware": self.sequence is not None and not self.legacy_library,
                "legacy_library": self.legacy_library, "expresses_alpha_helix": True,
                "expresses_beta_strand": True, "expresses_turns": True, "chiral": True,
                "realistic_bond_geometry": True}

LATTICE_DIRECTIONS = {(0, 0): (1.0, 1.0, 1.0), (0, 1): (1.0, -1.0, -1.0),
                      (1, 0): (-1.0, 1.0, -1.0), (1, 1): (-1.0, -1.0, 1.0)}
CA_CA_DISTANCE = 3.80
#: Lattice bond vectors have norm sqrt(3), so this puts adjacent CA atoms at 3.80 A.
LATTICE_SCALE = CA_CA_DISTANCE / math.sqrt(3.0)
_LATTICE_STEPS = np.array([LATTICE_DIRECTIONS[(a, b)]
                           for a in (0, 1) for b in (0, 1)], dtype=float)


class TetrahedralLatticeRepresentation:
    """CA-only tetrahedral lattice trace, 2 bits per bond."""

    name = "lattice"
    is_lattice = True

    def __init__(self, n_residues: int, sequence: Optional[str] = None):
        self.n_residues = int(n_residues)
        self.sequence = sequence.strip().upper() if sequence else None
        self.n_bonds = self.n_residues - 1
        self.n_bits = 2 * self.n_bonds
        self.n_qubits = self.n_bits
        self.bits_per_residue = 2

    def _check(self, bitstring: str) -> None:
        if len(bitstring) != self.n_bits:
            raise ValueError(f"bitstring length {len(bitstring)} != n_bits {self.n_bits}")

    def bond_directions(self, bitstring: str) -> List[Tuple[int, int]]:
        self._check(bitstring)
        return [(int(bitstring[2 * i]), int(bitstring[2 * i + 1]))
                for i in range(self.n_bonds)]

    def decode(self, bitstring: str) -> np.ndarray:
        """bitstring -> ``(n_res, 3)`` CA coordinates in ANGSTROMS.

        The `LATTICE_SCALE` factor was missing, and it was the largest single defect in the
        model: the raw lattice bond norm is sqrt(3) ~ 1.732, so every threshold in
        `energy_terms` -- all in Angstroms -- was wrong by 2.194x. The steric term charged
        20.80 for a |i-j| = 2 turn against 0.45 for a straight continuation, a 46x penalty
        for turning that came entirely from the unit mismatch and drove every search to the
        extended chain. Term-by-term accounting of the resulting all-straight optimum
        reproduces the previously published lattice minimum of 9.847041.
        """
        self._check(bitstring)
        d = np.frombuffer(bitstring.encode("latin-1"), dtype=np.uint8) - ord("0")
        idx = d[0::2] * 2 + d[1::2]
        out = np.zeros((self.n_residues, 3))
        np.cumsum(_LATTICE_STEPS[idx], axis=0, out=out[1:])
        return out * LATTICE_SCALE

    def build_coords(self, bitstring: str) -> Dict[str, np.ndarray]:
        """CA only. It used to also return ``"CB": ca.copy()``, which made every CB distance
        a duplicate of a CA distance and made the steric term charge the same pair twice at
        two different thresholds."""
        return {"CA": self.decode(bitstring)}

    def build_all(self, bitstring: str):
        """``(phi, psi, coords, rings)``. A CA-only trace has no torsions and no sidechains,
        so three of the four are None -- which is what makes `energy_terms` report torsion,
        both hbond terms and both cooperativity terms as exactly zero here."""
        return None, None, {"CA": self.decode(bitstring)}, None

    def native_bitstring(self, native_phi=None, native_psi=None, native_ca=None,
                         rng=None, iterations: int = 20000) -> str:
        """Best lattice trace for a native CA chain, by annealed search on CA-RMSD."""
        if native_ca is None:
            raise ValueError("lattice native_bitstring requires native_ca")
        if rng is None:
            rng = np.random.default_rng(0)
        native_ca = np.asarray(native_ca, dtype=float)

        def score(bits_list) -> float:
            # allow_scale stays off deliberately now that decode() returns Angstroms: a
            # free similarity scale flatters the lattice and makes its RMSD incomparable
            # with the torsion arm's.
            return geo.ca_rmsd(self.decode("".join(format(d, "02b") for d in bits_list)),
                               native_ca, allow_scale=False)

        cur = list(rng.integers(0, 4, size=self.n_bonds))
        cur_s = score(cur)
        best, best_s = list(cur), cur_s
        T0, T1 = 2.0, 1e-3
        for k in range(iterations):
            frac = k / max(1, iterations - 1)
            T = T0 * (1 - frac) + T1 * frac
            cand = list(cur)
            cand[int(rng.integers(0, self.n_bonds))] = int(rng.integers(0, 4))
            cs = score(cand)
            if cs < cur_s or rng.random() < math.exp(-(cs - cur_s) / max(T, 1e-9)):
                cur, cur_s = cand, cs
                if cs < best_s:
                    best, best_s = list(cand), cs
        return "".join(format(d, "02b") for d in best)

    def random_bitstring(self, rng: np.random.Generator) -> str:
        return "".join(format(int(d), "02b") for d in rng.integers(0, 4, size=self.n_bonds))

    def describe(self) -> Dict[str, object]:
        return {"name": self.name, "n_states": 4, "bits_per_residue": 2,
                "n_qubits": self.n_qubits, "config_space": 4.0 ** self.n_bonds,
                "ca_ca_distance": CA_CA_DISTANCE, "expresses_alpha_helix": False,
                "expresses_beta_strand": "approximately", "expresses_turns": "coarsely",
                "chiral": False, "realistic_bond_geometry": False}


def make_representation(kind: str, n_residues: int, n_states: int = 4,
                        sequence: Optional[str] = None, legacy_library: bool = False,
                        chi_bits: bool = True):
    """Build a representation.

    `sequence` is optional for backward compatibility but should ALWAYS be supplied for the
    torsion representation: without it every residue falls back to the GENERAL library,
    glycine/proline/pre-proline lose their class-specific states, and no chi1 bits are
    allocated. `legacy_library=True` restores the pre-2026-07-28 encoding entirely.
    """
    kind = kind.lower()
    if kind == "torsion":
        return TorsionStateRepresentation(n_residues, n_states=n_states, sequence=sequence,
                                          legacy_library=legacy_library, chi_bits=chi_bits)
    if kind == "lattice":
        return TetrahedralLatticeRepresentation(n_residues, sequence=sequence)
    raise ValueError(f"unknown representation {kind!r}")
