"""LANE L -- one pass over `prots/`, cached, so every `long40` candidate bank is cheap.

`s8/generate_univ` caches the ENTIRE length-n window universe per target.  That is
affordable at n = 13, where the universe is ~13,000 windows; at n = 45 it is 1.7 M windows
(~900 MB of coordinates per target) and it is not.  The window universe at long length is
therefore never materialised.  What is cached instead is the CORPUS -- each protein's
sequence, CA trace and torsions, parsed once -- and windows are cut from it on demand.

    python -m s32.s32_L_corpus build
"""
from __future__ import annotations

import glob
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(BASE, "s32", "results", "L_corpus.npz")


def build(verbose=True):
    from core import geometry as geo
    paths = sorted(glob.glob(os.path.join(BASE, "prots", "*.pdb")))
    pdbs, seqs, cas, phis, psis = [], [], [], [], []
    t0 = time.time()
    for k, p in enumerate(paths):
        try:
            seq, coords, phi, psi = geo.native_coords_from_pdb(p)
        except Exception:
            continue
        ca = np.asarray(coords["CA"], np.float32)
        if len(seq) < 9 or len(ca) != len(seq):
            continue
        pdbs.append(os.path.basename(p)[:-4].upper())
        seqs.append(seq)
        cas.append(ca)
        phis.append(np.asarray(phi, np.float32))
        psis.append(np.asarray(psi, np.float32))
        if verbose and (k + 1) % 3000 == 0:
            print(f"  {k+1}/{len(paths)} ({time.time()-t0:.0f}s)", flush=True)
    tmp = CACHE + f".{os.getpid()}.tmp.npz"        # unique temp: a shared one is not atomic
    np.savez(tmp, pdb=np.array(pdbs, dtype=object), seq=np.array(seqs, dtype=object),
             ca=np.array(cas, dtype=object), phi=np.array(phis, dtype=object),
             psi=np.array(psis, dtype=object))
    os.replace(tmp, CACHE)
    if verbose:
        print(f"{len(pdbs)} chains, {sum(len(s) for s in seqs)} residues -> {CACHE} "
              f"({time.time()-t0:.0f}s)", flush=True)
    return len(pdbs)


_C = {}


def load():
    """The corpus, plus the three per-chain masks every bank build needs.

    `code`, `ok` and `fin` are length-INDEPENDENT properties of a chain, so they are
    computed once here rather than once per (target, chain) pair.  At 45 targets x 13,751
    chains that is the difference between a minute and an afternoon.
    """
    if not _C:
        import s7.audit as audit
        from core import data as cdata
        AL = cdata.ALPHABET
        # The repository has TWO 20-letter orders that agree on three letters, and a
        # sprint was once voided by using the wrong one.  This asserts, rather than
        # assumes, that the canonical order is the one `audit.encode` and `audit.B62`
        # are indexed by -- the encoding these codes are about to be gathered with.
        assert np.array_equal(audit.encode(AL), np.arange(20)), "alphabet mismatch"
        z = np.load(CACHE, allow_pickle=True)
        _C["pdb"] = list(z["pdb"]); _C["seq"] = list(z["seq"])
        _C["ca"] = list(z["ca"]); _C["phi"] = list(z["phi"]); _C["psi"] = list(z["psi"])
        code, ok, fin = [], [], []
        for seq, ca, phi, psi in zip(_C["seq"], _C["ca"], _C["phi"], _C["psi"]):
            c = np.asarray(ca, float)
            step = np.linalg.norm(np.diff(c, axis=0), axis=1)
            ok.append((step > 3.5) & (step < 4.1))
            ph = np.asarray(phi, float); ps = np.asarray(psi, float)
            # A residue is usable if its code is a real amino acid and its interior
            # torsions are finite.  The two chain-terminal torsions are supplied by the
            # builder, so they are not required to be finite here.
            good = np.array([a in AL for a in seq], bool)
            f = np.isfinite(ph) | ~good
            g = np.isfinite(ps) | ~good
            f[0] = True; g[-1] = True
            fin.append(good & f & g)
            code.append(np.asarray([AL.index(a) if a in AL else 0
                                    for a in seq], np.int64))
        _C["code"] = code; _C["ok"] = ok; _C["fin"] = fin
    return _C


if __name__ == "__main__":
    build()
