"""Compact float16 cache of the FULL 1280-d ESM-2 per-residue embeddings.

Every sequence feature in this project passes through `esm_features.embed`, which projects
ESM-2's 1280 dimensions onto **32** principal components (`esm_pca.npz`: W is 1280x32).
That is a 40x compression sitting upstream of every structural prediction the project makes,
and it has never been ablated. Section 7 of s5/FINDINGS.md now gives a concrete accuracy
target -- distance MAE below ~1.5 A, from 2.22 A -- so it is worth knowing whether the
bottleneck is the model or its input.

Storing the raw embeddings for the 6,790 training sequences costs
95,039 residues x 1280 x 2 bytes = 243 MB in float16, against 1.5 GB for `esm_cache.npz`,
and lets any projection be tried without reloading it. float16 is safe here: these are
post-LayerNorm activations in a narrow range, and the downstream model normalises anyway.
"""
import os
import sys

sys.path.insert(0, os.path.abspath("."))
import numpy as np

OUT = "s5/esmraw.npz"


def build(force: bool = False) -> str:
    if os.path.exists(OUT) and not force:
        return OUT
    import esm_features as ef
    import fragment_db as fdb
    import peptide_db as db
    seqs = sorted({p.seq for p in db.load()} | {f.seq for f in fdb.load()})
    print(f"caching raw ESM for {len(seqs)} sequences", flush=True)
    keys, vals, offs = [], [], [0]
    for i, s in enumerate(seqs):
        e = ef.raw(s)[0].astype(np.float16)
        keys.append(s); vals.append(e); offs.append(offs[-1] + len(e))
        if (i + 1) % 1000 == 0:
            print(f"  {i+1}/{len(seqs)}", flush=True)
    flat = np.concatenate(vals, 0)
    tmp = OUT + ".tmp.npz"
    np.savez(tmp, keys=np.array(keys, dtype=object), flat=flat,
             offs=np.array(offs, dtype=np.int64))
    os.replace(tmp, OUT)
    print(f"wrote {OUT}: {flat.shape} float16, {flat.nbytes/1e6:.0f} MB", flush=True)
    return OUT


_C = None


def load():
    """dict sequence -> (n, 1280) float32."""
    global _C
    if _C is None:
        z = np.load(OUT, allow_pickle=True)
        k, f, o = z["keys"], z["flat"], z["offs"]
        _C = {str(s): f[o[i]:o[i + 1]].astype(np.float32) for i, s in enumerate(k)}
    return _C


if __name__ == "__main__":
    build(force="--force" in sys.argv)
