"""One-off: compact per-residue ESM-2 bank from the 1.5 GB esm_cache.npz.

Writes s12/cache/esm_bank.npz with, per sequence in the bank:
    pca32 (n, 32) float16  -- the shipped whitened PCA-32 (esm_pca.npz), per residue
    pca128 (n, 128) float16 -- a fresh 128-d PCA fitted here on all residues of the bank
    con (n, n) float16      -- ESM-2 contact probabilities
Run ONCE, alone (needs ~2.5 GB for the duration).  Everything else reads the compact file.
"""
import os, sys, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s12", "cache", "esm_bank.npz")

def build():
    z = np.load(os.path.join(ROOT, "esm_cache.npz"), allow_pickle=True)
    seqs = [str(s) for s in z["seqs"]]; vals = z["vals"]
    print("bank sequences", len(seqs), flush=True)
    p = np.load(os.path.join(ROOT, "esm_pca.npz")); mu, W, sc = p["mu"], p["W"], p["scale"]
    # fresh 128-d PCA on a residue subsample (memory-bounded)
    rng = np.random.default_rng(0)
    rows = []
    for k in rng.permutation(len(seqs))[:6000]:
        rows.append(np.asarray(vals[k][0], np.float32))
    X = np.vstack(rows); mu2 = X.mean(0)
    U, S, Vt = np.linalg.svd(X - mu2, full_matrices=False)
    W2 = Vt[:128].T; sc2 = S[:128] / np.sqrt(len(X))
    del X, rows, U, S, Vt
    pca32, pca128, con = [], [], []
    for k in range(len(seqs)):
        e = np.asarray(vals[k][0], np.float32); c = np.asarray(vals[k][1], np.float32)
        pca32.append((((e - mu) @ W) / np.maximum(sc, 1e-6)).astype(np.float16))
        pca128.append((((e - mu2) @ W2) / np.maximum(sc2, 1e-6)).astype(np.float16))
        con.append(c.astype(np.float16))
    np.savez_compressed(OUT, seqs=np.array(seqs, dtype=object),
                        pca32=np.array(pca32, dtype=object), pca128=np.array(pca128, dtype=object),
                        con=np.array(con, dtype=object), mu128=mu2, W128=W2, scale128=sc2)
    print("wrote", OUT, flush=True)

_BANK = None
def load():
    """dict seq -> (pca32 (n,32), pca128 (n,128), con (n,n)) as float32."""
    global _BANK
    if _BANK is None:
        z = np.load(OUT, allow_pickle=True)
        _BANK = {str(s): (np.asarray(a, np.float32), np.asarray(b, np.float32), np.asarray(c, np.float32))
                 for s, a, b, c in zip(z["seqs"], z["pca32"], z["pca128"], z["con"])}
    return _BANK

if __name__ == "__main__":
    build()
