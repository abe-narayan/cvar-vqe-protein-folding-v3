"""Structure database of short peptides, with leakage-controlled train/test separation.

Why this exists. Every learned component in this project -- the torsion library, the
per-residue basin prior, the pairwise distance prior, the ranking model -- was previously
fitted on the 13 peptides in `dataset.CANDIDATE_PDB_IDS` plus whatever happened to be in
`pdbs/`. Thirteen targets is not a training set; `learned_fold`'s own docstring records
that pooling 45 more peptides made the scorer *worse*, which is the signature of a model
with too little data to separate fold-family from generality rather than evidence that
generality is impossible.

`pdbs_ext/` holds every RCSB entry that is a single protein chain of 8-26 residues
(1,470 entries, queried 2026-08-31). After parsing, torsion-rebuild fidelity gating and
exact-sequence deduplication, 787 usable structures remain -- a 60x increase.

Leakage control. `holdout` returns the database with every entry whose global alignment
identity to the target exceeds `IDENTITY_THRESHOLD` removed, plus the target itself. The
identity is the same Needleman-Wunsch measure `dataset._identity` uses, normalised by the
longer sequence, so a short peptide fully contained in a longer one does not score 1.00
and slip through. Identity is computed target-vs-database only (a few thousand pairs), not
all-vs-all, because that is the comparison leakage actually depends on.
"""
import glob
import json
import os
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import protein_geometry as geo

BASE = os.path.dirname(os.path.abspath(__file__))
DIRS = (os.path.join(BASE, "pdbs"), os.path.join(BASE, "pdbs_ext"))
CACHE = os.path.join(BASE, "peptide_db.npz")

MIN_LEN, MAX_LEN = 8, 26
#: Backbone (N, CA, C) RMSD between a native and the structure rebuilt from its own
#: torsions by the ideal-geometry builder. Above this the deposited geometry carries
#: something the representation cannot express (cis peptide bonds, distorted valence),
#: so the entry is not evidence about anything the search can reach.
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
            pdb = os.path.basename(path)[:-4].upper()
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
                seen[seq] = dict(pdb=pdb, seq=seq, ca=ca, phi=phi, psi=psi, rebuild=rb)
    return sorted(seen.values(), key=lambda r: r["pdb"])


def build(force: bool = False) -> List[Peptide]:
    if not force and os.path.exists(CACHE):
        z = np.load(CACHE, allow_pickle=True)
        return [Peptide(**r) for r in z["records"]]
    recs = _scan()
    out = [Peptide(**r) for r in recs]
    np.savez_compressed(CACHE, records=np.array(recs, dtype=object))
    return out


@lru_cache(maxsize=1)
def load() -> Tuple[Peptide, ...]:
    return tuple(build())


def by_pdb(pdb: str) -> Optional[Peptide]:
    pdb = pdb.upper()
    for p in load():
        if p.pdb == pdb:
            return p
    return None


def identity(a: str, b: str, gap: float = -1.0) -> float:
    """Needleman-Wunsch identity normalised by the longer sequence (see `dataset`)."""
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return 0.0
    prev = [(gap * j, 0) for j in range(m + 1)]
    for i in range(1, n + 1):
        cur = [(gap * i, 0)]
        ai = a[i - 1]
        for j in range(1, m + 1):
            match = 1.0 if ai == b[j - 1] else 0.0
            diag = (prev[j - 1][0] + match, prev[j - 1][1] + int(match > 0))
            up = (prev[j][0] + gap, prev[j][1])
            left = (cur[j - 1][0] + gap, cur[j - 1][1])
            best = diag
            if up[0] > best[0]:
                best = up
            if left[0] > best[0]:
                best = left
            cur.append(best)
        prev = cur
    return prev[m][1] / max(n, m)


def _kmers(s: str, k: int = 3) -> set:
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def _composition(s: str) -> Dict[str, int]:
    c: Dict[str, int] = {}
    for ch in s:
        c[ch] = c.get(ch, 0) + 1
    return c


def max_possible_identity(a: str, b: str,
                          ca: Optional[Dict[str, int]] = None,
                          cb: Optional[Dict[str, int]] = None) -> float:
    """An ADMISSIBLE upper bound on `identity(a, b)`.

    No alignment can match a residue type more often than it occurs in the scarcer of the
    two sequences, so the multiset intersection bounds the match count from above, and
    `identity` divides by the longer length. This is a sound prefilter -- it can only
    over-estimate -- unlike the 3-mer Jaccard test it replaces, whose premise ("identity
    above 0.6 forces a shared 3-mer") is false because matches may alternate with
    mismatches. That unsound filter let one benchmark target train on a 0.70-identity
    homolog. This bound is 20 dictionary lookups and skips ~97% of pairs.
    """
    ca = ca or _composition(a)
    cb = cb or _composition(b)
    shared = sum(min(v, cb.get(ch, 0)) for ch, v in ca.items())
    return shared / max(len(a), len(b))


def holdout(target_seq: str, threshold: float = IDENTITY_THRESHOLD
            ) -> Tuple[Peptide, ...]:
    """Database with every entry too similar to `target_seq` removed.

    EXACT: every database entry is aligned against the target. The 3-mer prefilter this
    used to carry was unsound -- see `clusters` -- and one benchmark target was provably
    trained on a 0.70-identity homolog because of it. 787 alignments per target is cheap.
    """
    target = target_seq.strip().upper()
    ct = _composition(target)
    out = []
    for p in load():
        if p.seq == target:
            continue
        if max_possible_identity(target, p.seq, ct) < threshold:
            out.append(p)
            continue
        if identity(target, p.seq) < threshold:
            out.append(p)
    return tuple(out)


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
    # Deduplication may have replaced a local entry with a better-rebuilding deposit of
    # the SAME sequence (1J4M -> 1K43, 5AWL -> 2RVD). Selecting by sequence keeps the
    # target rather than silently dropping it.
    return [p for p in load()
            if p.seq in keep_seq and min_len <= p.n <= max_len]


if __name__ == "__main__":
    db = build(force=True)
    print(f"{len(db)} peptides, lengths {min(p.n for p in db)}-{max(p.n for p in db)}")
    print(f"median rebuild RMSD {np.median([p.rebuild for p in db]):.2f} A")
    bench = benchmark_set()
    print(f"benchmark: {len(bench)} targets -> "
          f"{', '.join(sorted(p.pdb for p in bench))}")


# --------------------------------------------------------------- clustering / folds
_CLUSTER_CACHE = os.path.join(BASE, "peptide_clusters.json")


def clusters(threshold: float = IDENTITY_THRESHOLD) -> Dict[str, int]:
    """Single-linkage identity clusters over the whole database, sequence -> cluster id.

    EXACT: every pair is aligned. The previous version skipped pairs with no shared 3-mer,
    justified by "a pair above 0.6 identity must share one" -- which is false. Identity is
    a count of matches along an alignment and those matches can alternate with mismatches,
    so a pair can reach 0.6 with no shared 3-mer at all. A brute-force pass over all
    309,141 pairs found two such pairs, one of them real and consequential: 1C9A, a
    reported benchmark target, was trained on 2N0G at 0.70 identity. The saving was never
    worth an unsound guarantee; the exact pass takes a few minutes and is cached.
    """
    if os.path.exists(_CLUSTER_CACHE):
        with open(_CLUSTER_CACHE) as f:
            saved = json.load(f)
        if saved.get("threshold") == threshold:
            return saved["assign"]
    seqs = [p.seq for p in load()]
    comps = [_composition(s) for s in seqs]
    parent = list(range(len(seqs)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a in range(len(seqs)):
        for b in range(a + 1, len(seqs)):
            if find(a) == find(b):
                continue
            if max_possible_identity(seqs[a], seqs[b], comps[a], comps[b]) < threshold:
                continue
            if identity(seqs[a], seqs[b]) >= threshold:
                parent[find(b)] = find(a)
    roots = {}
    assign = {}
    for i, s in enumerate(seqs):
        r = find(i)
        assign[s] = roots.setdefault(r, len(roots))
    with open(_CLUSTER_CACHE, "w") as f:
        json.dump({"threshold": threshold, "assign": assign}, f)
    return assign


_FOLD_CACHE = os.path.join(BASE, "peptide_folds.json")


def folds(n_folds: int = 5, seed: int = 0) -> Dict[str, int]:
    """Sequence -> fold index, assigning whole identity clusters to folds. PINNED.

    Cross-validating over folds of clusters rather than per target costs one model per
    fold instead of one per target, with the same guarantee: a target's fold model never
    saw the target or anything above the identity threshold to it.

    The assignment is written to `peptide_folds.json` on first use and read back
    thereafter, and that pinning is load-bearing rather than a cache. Fold indices are a
    shuffle of CLUSTER IDS, so anything that renumbers clusters -- adding a database entry,
    or replacing an unsound identity prefilter with a sound one -- silently reassigns every
    sequence. The trained fold models on disk would then be excluding the wrong fold, and a
    target could be scored by a model that had trained on it, with no error raised
    anywhere. Delete this file only together with every model in `distogram_models/`,
    `pairnet_models/` and `torsion_models/`.
    """
    if os.path.exists(_FOLD_CACHE):
        with open(_FOLD_CACHE) as f:
            saved = json.load(f)
        if saved.get("n_folds") == n_folds and saved.get("seed") == seed:
            assign = saved["assign"]
            missing = [p.seq for p in load() if p.seq not in assign]
            if not missing:
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


def dev_set(n: int = 24, min_len: int = 9, max_len: int = 16,
            seed: int = 7, _seq_only: bool = False) -> List[Peptide]:
    """Hyperparameter development targets, disjoint from `benchmark_set`.

    Objective weights, CVaR schedule constants and selection radii are chosen on these and
    never on the reported targets. The database is large enough that this costs nothing:
    752 peptides are outside the benchmark. One peptide per identity cluster, so the dev
    set is not 24 views of the same fold.
    """
    # CLUSTER-disjoint, not merely sequence-disjoint. The earlier version excluded only
    # exact benchmark sequences, which left four development peptides sharing an identity
    # cluster with six reported targets (1Y58 / 1Y5C at 0.846 identity; 1B45 / 1R8T and
    # 1QMW; 1AKG / 1UL2, 2H8S, 2JUT). Hyperparameters fitted on that set measurably
    # transferred to its homologs and not to the rest: a controlled re-run of the fitted
    # score weights helped the six contaminated targets by 0.135 A and hurt eight controls
    # by 0.289 A, a difference-in-differences of 0.424 A. Anything fitted here is fitted on
    # clusters the benchmark never contains.
    # `_seq_only=True` is the internal, acyclic form: `benchmark()` needs a dev set to
    # exclude while it is being defined, so it asks for the sequence-only version. The
    # public form then excludes the whole benchmark's CLUSTERS, which is what makes
    # anything fitted here safe to report on the benchmark. Benchmark membership is fixed
    # by the internal form and does not depend on this stronger filter, so it is stable.
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


#: The reported benchmark. Membership is fixed by this function and by
#: `results/benchmark_manifest.json`; it is never filtered on difficulty or on any result.
#:
#: Construction, in order, so it is reproducible and auditable:
#:   1. every peptide in the repository's own `pdbs/` directory in the length range --
#:      the 35 targets the earlier work reported, kept so the comparison is continuous;
#:   2. cluster representatives drawn from the rest of the database by a seeded shuffle
#:      over identity clusters, one per cluster, until `n` targets are reached.
#: The `dev_set` peptides are excluded at both steps: they carry the hyperparameters.
BENCHMARK_N = 60
BENCHMARK_MIN_LEN = 9
BENCHMARK_MAX_LEN = 16


MANIFEST = os.path.join(BASE, "results", "benchmark_manifest.json")
#: The MONOMER benchmark: 52 peptides whose deposited conformation is that of an isolated
#: monomer. The 60-target set was shown to be dominated by chains that are not folded on
#: their own -- fibril segments (7.22 A), micelle- and membrane-bound peptides (5.27 A),
#: partner-bound (4.03 A) -- where the reference conformation is held by inter-chain packing
#: or by a detergent, and a single-chain intramolecular energy has no access to it. The
#: monomer set is selected on the EXPERIMENT (deposited header text, model count, ensemble
#: spread, builder reproducibility, cluster uniqueness), never on any method's performance,
#: and `work/monomer_benchmark.md` logs every one of the 256 exclusions with the deposited
#: line that triggered it. 19 of the 60 survive into it.
MONOMER_MANIFEST = os.path.join(BASE, "results", "monomer_manifest.json")


def monomer_benchmark(aqueous_only: bool = False) -> List[Peptide]:
    """The 52-target monomer benchmark, or the 36-target strictly-aqueous subset.

    `aqueous_only=True` additionally drops peptides solved in 30-100% TFE, HFIP, DMSO or
    methanol: those are genuine isolated monomers, so they meet the criterion as written,
    but their conformation is solvent-imposed. The subset is nested inside the 52.
    """
    with open(MONOMER_MANIFEST) as f:
        man = json.load(f)
    by_id = {p.pdb: p for p in load()}
    out = []
    for t in man["targets"]:
        if aqueous_only and not t.get("aqueous", True):
            continue
        if t["pdb"] in by_id:
            out.append(by_id[t["pdb"]])
    return out


def benchmark(n: int = BENCHMARK_N, min_len: int = BENCHMARK_MIN_LEN,
              max_len: int = BENCHMARK_MAX_LEN, seed: int = 11) -> List[Peptide]:
    """The reported benchmark, FROZEN to `results/benchmark_manifest.json` once written.

    Membership is generated once (see below) and then pinned. This is not a convenience:
    the selection below depends on identity-cluster ids, and a change to the clustering
    -- for example replacing an unsound k-mer prefilter with a sound one, which merged one
    pair and renumbered everything downstream -- silently swapped 13 of the 60 targets.
    A benchmark whose membership moves when an internal detail changes cannot support a
    before/after claim, and re-running the baseline is not the fix; freezing is.

    Delete the manifest to regenerate deliberately.
    """
    if os.path.exists(MANIFEST):
        with open(MANIFEST) as f:
            want = [t["pdb"] for t in json.load(f)["targets"]]
        by_id = {p.pdb: p for p in load()}
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
