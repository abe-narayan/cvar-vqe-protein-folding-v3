"""Build a compact 32-d ESM cache for every sequence Sprint 5 trains on.

`esm_features.raw` loads `esm_cache.npz`, which is 1.5 GB, into a module global. With a
few hundred MB free on this box that is a swap event, and it happens on the first call in
every process that touches ESM. The projected embedding is only 32 dimensions, so the
whole training set fits in about 13 MB:

    6,800 structures x ~14 residues x 32 float32

This builds that once, in a short-lived process that exits and returns the 1.5 GB, and
every later process memory-maps the small file instead.
"""
import os
import sys

sys.path.insert(0, os.path.abspath("."))
import numpy as np

OUT = "s5/esm32.npz"


def build(force: bool = False) -> str:
    if os.path.exists(OUT) and not force:
        return OUT
    import distogram as dgm
    import esm_features as ef
    import fragment_db as fdb
    import peptide_db as db
    seqs = sorted({p.seq for p in db.load()} | {f.seq for f in fdb.load()})
    print(f"embedding {len(seqs)} unique sequences", flush=True)
    keys, vals, offs = [], [], [0]
    for i, s in enumerate(seqs):
        e = ef.embed(s).astype(np.float32)
        keys.append(s)
        vals.append(e)
        offs.append(offs[-1] + len(e))
        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{len(seqs)}", flush=True)
    flat = np.concatenate(vals, 0)
    tmp = OUT + ".tmp.npz"
    np.savez_compressed(tmp, keys=np.array(keys, dtype=object),
                        flat=flat, offs=np.array(offs, dtype=np.int64))
    os.replace(tmp, OUT)
    print(f"wrote {OUT}: {flat.shape} ({flat.nbytes/1e6:.1f} MB raw)", flush=True)
    return OUT


_CACHE = None


def load():
    """dict sequence -> (n, 32) float32. ~13 MB resident."""
    global _CACHE
    if _CACHE is None:
        z = np.load(OUT, allow_pickle=True)
        keys, flat, offs = z["keys"], z["flat"], z["offs"]
        _CACHE = {str(k): flat[offs[i]:offs[i + 1]] for i, k in enumerate(keys)}
    return _CACHE


if __name__ == "__main__":
    build(force="--force" in sys.argv)
