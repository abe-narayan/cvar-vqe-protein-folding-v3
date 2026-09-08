"""One-off, run ALONE: extract the 787 peptide-database sequences from the compact
`s12/cache/esm_bank.npz` into a tiny per-peptide bank so that every later process can
serve `core.data.esm_embed` / `esm_contacts` without touching either the 1.5 GB
`esm_cache.npz` (forbidden) or the 100 MB bank.

Writes s12/cache/dir_esm_pep.npz  (seqs, pca32 float32, con float32) ~ 6 MB.
"""
import os, sys, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s12", "cache", "dir_esm_pep.npz")


def main():
    recs = np.load(os.path.join(ROOT, "peptide_db.npz"), allow_pickle=True)["records"]
    want = list(dict.fromkeys(r["seq"] for r in recs))
    z = np.load(os.path.join(ROOT, "s12", "cache", "esm_bank.npz"), allow_pickle=True)
    seqs = [str(s) for s in z["seqs"]]
    idx = {s: k for k, s in enumerate(seqs)}
    pca32, con, keep = [], [], []
    P, C = z["pca32"], z["con"]
    for s in want:
        if s not in idx:
            continue
        k = idx[s]
        keep.append(s)
        pca32.append(np.asarray(P[k], np.float32))
        con.append(np.asarray(C[k], np.float32))
    np.savez_compressed(OUT, seqs=np.array(keep, dtype=object),
                        pca32=np.array(pca32, dtype=object), con=np.array(con, dtype=object))
    print("wrote", OUT, len(keep), "of", len(want), flush=True)


if __name__ == "__main__":
    main()
