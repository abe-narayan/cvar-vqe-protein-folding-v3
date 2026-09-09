"""ESM-2 features for the peptide database, cached to disk.

What is taken from the model and why. Two things, both per sequence:

* **Per-residue representations** from the final layer, reduced to `N_PCA` principal
  components fitted once over the whole database. A 1280-dimensional embedding against
  ~13,000 database residues would let the downstream tree model memorise identity rather
  than learn a mapping; the PCA is a fixed, sequence-independent projection that keeps the
  directions the peptide corpus actually varies along.
* **Attention-derived contact logits** from ESM-2's own contact head, which is a logistic
  regression over symmetrised, APC-corrected attention maps trained on real contacts. That
  is a pairwise quantity, which is exactly the shape the distance prior needs, and it is
  the only part of ESM-2 that was ever supervised on structure.

Both are cached in `esm_cache.npz` keyed by sequence, so a benchmark sweep pays the model
cost once. Nothing here reads a structure, so ESM features are available for a target whose
native is being held out.
"""
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "esm_cache.npz")
#: A small hot cache holding only the sequences inference actually needs. The full cache
#: grew to 1.5 GB once 22,000 fragments were featurised, and `np.load` on an object-array
#: npz materialises all of it -- which a 16 GB machine running a benchmark cannot afford.
#: Training reads `CACHE`; everything else reads this first and only falls back.
CACHE_SMALL = os.path.join(BASE, "esm_small.npz")
MODEL_NAME = "esm2_t33_650M_UR50D"
N_PCA = 32

_cache: Optional[Dict[str, Tuple[np.ndarray, np.ndarray]]] = None
_small: Optional[Dict[str, Tuple[np.ndarray, np.ndarray]]] = None


def _read(path: str) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    if not os.path.exists(path):
        return {}
    z = np.load(path, allow_pickle=True)
    return {str(k): (v[0], v[1]) for k, v in zip(z["seqs"], z["vals"])}


def _load_small() -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    global _small
    if _small is None:
        _small = _read(CACHE_SMALL)
    return _small


def _load_cache() -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    global _cache
    if _cache is None:
        _cache = _read(CACHE)
    return _cache


def write_small(sequences: Sequence[str]) -> None:
    """Featurise `sequences` (computing any that are missing) into the hot cache."""
    have = _load_small()
    todo = [s for s in dict.fromkeys(sequences) if s not in have]
    if todo:
        compute(todo, verbose=False)
        big = _load_cache()
        have.update({s: big[s] for s in todo if s in big})
    keys = list(have.keys())
    np.savez_compressed(CACHE_SMALL, seqs=np.array(keys, dtype=object),
                        vals=np.array([[have[s][0], have[s][1]] for s in keys],
                                      dtype=object))


def available() -> bool:
    return os.path.exists(CACHE) or os.path.exists(CACHE_SMALL)


def compute(sequences: Sequence[str], batch: int = 8, verbose: bool = True) -> None:
    """Run ESM-2 over `sequences` and extend the cache. Requires torch + fair-esm."""
    import torch
    import esm

    cache = _load_cache()
    todo = [s for s in dict.fromkeys(sequences) if s not in cache]
    if not todo:
        return
    model, alphabet = getattr(esm.pretrained, MODEL_NAME)()
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
                cache[s] = (reps[m, 1:n + 1].astype(np.float32),
                            cons[m, :n, :n].astype(np.float32))
            if verbose:
                print(f"  esm {start + len(chunk)}/{len(todo)}", flush=True)
    seqs = list(cache.keys())
    vals = np.array([[cache[s][0], cache[s][1]] for s in seqs], dtype=object)
    np.savez_compressed(CACHE, seqs=np.array(seqs, dtype=object), vals=vals)


def raw(sequence: str) -> Tuple[np.ndarray, np.ndarray]:
    """``(per-residue representation (n, D), contact probabilities (n, n))``."""
    hot = _load_small()
    if sequence in hot:
        return hot[sequence]
    cache = _load_cache()
    if sequence not in cache:
        compute([sequence], verbose=False)
        cache = _load_cache()
    return cache[sequence]

# ------------------------------------------------------------------ projection
_PCA: Optional[Tuple[np.ndarray, np.ndarray]] = None
_PCA_PATH = os.path.join(BASE, "esm_pca.npz")


def fit_pca(sequences: Sequence[str], n_components: int = N_PCA) -> None:
    X = np.vstack([raw(s)[0] for s in sequences])
    mu = X.mean(0)
    U, S, Vt = np.linalg.svd(X - mu, full_matrices=False)
    np.savez(_PCA_PATH, mu=mu, W=Vt[:n_components].T,
             scale=(S[:n_components] / np.sqrt(len(X))))


def _pca() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    global _PCA
    if _PCA is None:
        z = np.load(_PCA_PATH)
        _PCA = (z["mu"], z["W"], z["scale"])
    return _PCA


def embed(sequence: str) -> np.ndarray:
    """``(n, N_PCA)`` whitened per-residue embedding."""
    mu, W, scale = _pca()
    return ((raw(sequence)[0] - mu) @ W) / np.maximum(scale, 1e-6)


def contacts(sequence: str) -> np.ndarray:
    return raw(sequence)[1]
