"""Build the leave-fold-out training corpora for the torsion-bin predictor.

The universe file of a target T contains EXACTLY the legal library for T (out-of-fold
peptides + this fold's fragments) as every length-n window of every parent chain, in
parent order and sliding by one residue.  Consecutive windows therefore overlap by n-1,
which lets us STITCH the parents back and recover, per parent chain:

    codes (L,) int8 amino-acid codes,  phi (L,), psi (L,) radians,  org bool

This gives a redundancy-free per-residue training set that is leakage-safe by construction
for every target in that fold.  Cached to s12/cache/key_corpus_f<fold>.npz.
"""
from __future__ import annotations
import os, sys
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I


def stitch(u):
    """Recover parent chains from the sliding windows of a universe."""
    S = np.asarray(u["S"], np.int8); PH = np.asarray(u["PHI"], float); PS = np.asarray(u["PSI"], float)
    org = np.asarray(u["org"], bool)
    nw, n = S.shape
    cont = np.zeros(nw, bool)
    if nw > 1:
        cont[1:] = (S[1:, :-1] == S[:-1, 1:]).all(1) & (org[1:] == org[:-1]) & \
                   (np.abs(PH[1:, :-1] - PH[:-1, 1:]) < 1e-3).all(1)
    starts = np.where(~cont)[0]
    ends = np.append(starts[1:], nw)
    out = []
    for a, b in zip(starts, ends):
        L = n + (b - a - 1)
        codes = np.empty(L, np.int8); phi = np.empty(L); psi = np.empty(L)
        codes[:n] = S[a]; phi[:n] = PH[a]; psi[:n] = PS[a]
        if b - a > 1:
            codes[n:] = S[a + 1:b, -1]; phi[n:] = PH[a + 1:b, -1]; psi[n:] = PS[a + 1:b, -1]
        out.append((codes, phi, psi, bool(org[a])))
    return out


def _fold_source(fold):
    """The universe to stitch for a fold: the smallest-n target in that fold
    (smallest n recovers the most parents and the longest chains)."""
    tg = [t for t in I.targets() if t["fold"] == fold]
    return min(tg, key=lambda t: (t["n"], t["pdb"]))["pdb"]


def corpus(fold):
    """dict codes/phi/psi/org/pid arrays (concatenated parents) for `fold`'s legal library."""
    path = os.path.join(I.CACHE, f"key_corpus_f{fold}.npz")
    if os.path.exists(path):
        z = np.load(path)
        return {k: z[k] for k in z.files}
    pdb = _fold_source(fold)
    u = I.load_univ(pdb)
    ps = stitch(u)
    codes = np.concatenate([p[0] for p in ps])
    phi = np.concatenate([p[1] for p in ps]).astype(np.float32)
    psi = np.concatenate([p[2] for p in ps]).astype(np.float32)
    org = np.concatenate([np.full(len(p[0]), p[3]) for p in ps])
    pid = np.concatenate([np.full(len(p[0]), k, np.int32) for k, p in enumerate(ps)])
    out = {"codes": codes, "phi": phi, "psi": psi, "org": org, "pid": pid,
           "src": np.array([pdb])}
    np.savez_compressed(path, **out)
    del u
    return out


if __name__ == "__main__":
    for f in sorted({t["fold"] for t in I.targets()}):
        c = corpus(f)
        npar = len(np.unique(c["pid"]))
        print(f"fold {f}: src={c['src'][0]} parents={npar} residues={len(c['codes'])} "
              f"peptide-res={int(c['org'].sum())} frag-res={int((~c['org']).sum())} free={I.free_gb():.2f}",
              flush=True)
