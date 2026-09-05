"""Backbone fragments of small proteins, as extra training data for the distance prior.

The distance prior is the component that decides final accuracy (see `distogram`), and it
was being fitted on 640 peptides -- about 85,000 residue pairs -- for a model with 183
inputs. Contiguous windows of small, high-resolution crystal structures are the obvious
place to get more of the same quantity: the geometry of a nine-residue stretch of chain is
the same object whether or not the rest of a protein is attached to it.

The caveat is real and is why fragments are training-only. A fragment's conformation is
partly held in place by contacts outside the window, so the fragment distribution is more
compact and more regular than the distribution of *isolated* peptides. It is a source of
geometry, not a source of "what does this peptide do on its own", and the peptide database
remains what the torsion library, the MRF prior and the benchmark are built from.
`work/eval_distogram_frag.py` measures the effect rather than assuming it.

Leakage control is the same as everywhere else: `holdout` drops any fragment above the
identity threshold to the target.
"""
import glob
import os
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import peptide_db as pdb
import protein_geometry as geo

BASE = os.path.dirname(os.path.abspath(__file__))
PROT_DIR = os.path.join(BASE, "prots")
CACHE = os.path.join(BASE, "fragment_db.npz")
#: A second, larger extraction kept in its own file so the 6,003-fragment result
#: stays reproducible while a bigger set is evaluated against it.
CACHE_LARGE = os.path.join(BASE, "fragment_db_large.npz")

LENGTHS = tuple(range(9, 21))
STRIDE = 5
REBUILD_TOL = 1.0


def _extract(path: str) -> List[Dict]:
    try:
        seq, coords, phi, psi = geo.native_coords_from_pdb(path)
    except Exception:
        return []
    ca = np.asarray(coords["CA"], float)
    n = len(seq)
    if n < min(LENGTHS) or len(ca) != n:
        return []
    step = np.linalg.norm(np.diff(ca, axis=0), axis=1)
    ok = (step > 3.5) & (step < 4.1)
    out = []
    for L in LENGTHS:
        for s in range(1, n - L, STRIDE):        # skip residue 0: phi is undefined
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
            out.append(dict(pdb=os.path.basename(path)[:-4].upper() + f"_{s}",
                            seq=sub, ca=ca[s:e].copy(), phi=ph.copy(), psi=ps.copy(),
                            rebuild=0.0))
    return out


def build(max_per_protein: int = 6, n_target: int = 6000, seed: int = 0,
          force: bool = False, cache: Optional[str] = None) -> List[pdb.Peptide]:
    cache = cache or CACHE
    if not force and os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return [pdb.Peptide(**r) for r in z["records"]]
    rng = np.random.default_rng(seed)
    seen: Dict[str, Dict] = {}
    paths = sorted(glob.glob(os.path.join(PROT_DIR, "*.pdb")))
    rng.shuffle(paths)
    for path in paths:
        frags = _extract(path)
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
                m = geo.kabsch_superpose(built["CA"], f["ca"])
                if geo.rmsd(m, f["ca"]) > REBUILD_TOL:
                    continue
            except Exception:
                continue
            seen[f["seq"]] = f
            kept += 1
        if len(seen) >= n_target:
            break
    recs = list(seen.values())
    np.savez_compressed(cache, records=np.array(recs, dtype=object))
    return [pdb.Peptide(**r) for r in recs]


#: Set to True to make every consumer read the larger extraction. Kept as a switch rather
#: than a silent default so a result can always be attributed to a fragment set.
USE_LARGE = os.environ.get("FRAG_LARGE", "") == "1"


@lru_cache(maxsize=2)
def load(large: Optional[bool] = None) -> Tuple[pdb.Peptide, ...]:
    if large is None:
        large = USE_LARGE
    return tuple(build(cache=CACHE_LARGE if large else CACHE))


def holdout(target_seq: str,
            threshold: float = pdb.IDENTITY_THRESHOLD) -> Tuple[pdb.Peptide, ...]:
    target = target_seq.strip().upper()
    tk = pdb._kmers(target)
    out = []
    for p in load():
        if p.seq == target:
            continue
        if tk and not (pdb._kmers(p.seq) & tk):
            out.append(p)
            continue
        if pdb.identity(target, p.seq) >= threshold:
            continue
        out.append(p)
    return tuple(out)


if __name__ == "__main__":
    frags = build(force=True)
    import collections
    print(f"{len(frags)} fragments from {len(set(f.pdb.split('_')[0] for f in frags))} "
          f"proteins")
    print(collections.Counter(f.n for f in frags))
