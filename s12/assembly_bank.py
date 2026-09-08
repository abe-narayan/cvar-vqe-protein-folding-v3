"""Sprint 12 / assembly: per-(fold, L) short-piece banks from the leakage-safe library.

For fold f the library is exactly what the pipeline's universe uses: out-of-fold peptides
(peptide_db, folds(5)[seq] != f) plus this fold's fragment set (core.predict._fold_fragments(f, 5)).

Residue table  s12/cache/pieces_f{f}_res.npz
    phi, psi (R,) float32 radians     code (R,) int8      chain (R,) int32   pos (R,) int16
    clen (R,) int16   is_pep (R,) bool   esm (R,128) float16 (ESM-2 pca128 from s12/esm_bank)
    valid_phi/valid_psi (R,) bool  -- False where the value is a chain-end placeholder
    state (R,) int8  -- circular k-means (k=8) torsion state fitted on THIS fold's library residues
    centres (8,2)    -- state centres (phi, psi)
    trans (8,8)      -- bigram P(state_{i+1} | state_i) over consecutive residues (add-1 smoothed)
Piece bank     s12/cache/pieces_f{f}_L{L}.npz
    start (P,) int32  -- residue-table index of the first residue; piece = rows start..start+L-1
    ca (P,L,3) float32 -- ideal-geometry CA built from the piece's own torsions (instrument.build_ca)
    rca (P,L,3) float32 -- the parent's REAL CA for the same residues
"""
from __future__ import annotations
import os, sys, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I

LS = list(range(3, 11))
NSTATE = 8


def _chains(fold):
    import peptide_db as db
    from core import predict
    folds = db.folds(5)
    peps = [q for q in db.load() if folds[q.seq] != fold]
    frs = list(predict._fold_fragments(fold, 5))
    return peps, frs


def _circ_kmeans(X, k=NSTATE, seed=0, iters=50):
    """k-means on the unit-circle embedding of (phi, psi)."""
    E = np.concatenate([np.cos(X), np.sin(X)], 1)
    rng = np.random.default_rng(seed)
    # deterministic init from canonical regions so the state ids are stable across folds
    init = np.deg2rad(np.array([[-63, -42], [-120, 130], [-70, 140], [-90, 0], [60, 40], [-140, 160], [-75, -20], [80, -170]], float))
    C = np.concatenate([np.cos(init), np.sin(init)], 1)
    for _ in range(iters):
        d = ((E[:, None, :] - C[None]) ** 2).sum(-1)
        a = d.argmin(1)
        for c in range(k):
            m = a == c
            if m.any():
                C[c] = E[m].mean(0)
    C4 = C.reshape(k, 2, 2)
    centres = np.arctan2(C4[:, 1, :], C4[:, 0, :])
    return centres, a


def assign_state(phi, psi, centres):
    E = np.stack([np.cos(phi), np.cos(psi), np.sin(phi), np.sin(psi)], -1)
    C = np.concatenate([np.cos(centres), np.sin(centres)], 1)
    d = ((E[..., None, :] - C) ** 2).sum(-1)
    return d.argmin(-1)


def res_path(fold):
    return os.path.join(I.CACHE, f"pieces_f{fold}_res.npz")


def piece_path(fold, L):
    return os.path.join(I.CACHE, f"pieces_f{fold}_L{L}.npz")


def build_fold(fold, verbose=True):
    from s12 import esm_bank
    if all(os.path.exists(piece_path(fold, L)) for L in LS) and os.path.exists(res_path(fold)):
        return
    peps, frs = _chains(fold)
    bank = esm_bank.load()
    phi, psi, code, chain, pos, clen, is_pep, esm, vphi, vpsi = [], [], [], [], [], [], [], [], [], []
    ci = 0
    for kind, lst in (("pep", peps), ("frag", frs)):
        for q in lst:
            m = q.n
            e = bank.get(q.seq)
            if e is None:
                raise RuntimeError(f"{q.pdb} missing from esm bank")
            phi.append(np.asarray(q.phi, np.float32)); psi.append(np.asarray(q.psi, np.float32))
            code.append(np.array([I.ALPHABET.index(c) for c in q.seq], np.int8))
            chain.append(np.full(m, ci, np.int32)); pos.append(np.arange(m, dtype=np.int16))
            clen.append(np.full(m, m, np.int16)); is_pep.append(np.full(m, kind == "pep", bool))
            esm.append(np.asarray(e[1], np.float16))
            vp = np.ones(m, bool); vs = np.ones(m, bool)
            if kind == "pep":
                vp[0] = False; vs[-1] = False
            vphi.append(vp); vpsi.append(vs)
            ci += 1
    R = {"phi": np.concatenate(phi), "psi": np.concatenate(psi), "code": np.concatenate(code),
         "chain": np.concatenate(chain), "pos": np.concatenate(pos), "clen": np.concatenate(clen),
         "is_pep": np.concatenate(is_pep), "esm": np.concatenate(esm),
         "valid_phi": np.concatenate(vphi), "valid_psi": np.concatenate(vpsi)}
    ok = R["valid_phi"] & R["valid_psi"]
    centres, _ = _circ_kmeans(np.stack([R["phi"][ok], R["psi"][ok]], 1).astype(float))
    st = assign_state(R["phi"].astype(float), R["psi"].astype(float), centres).astype(np.int8)
    R["state"] = st; R["centres"] = centres
    T = np.ones((NSTATE, NSTATE))
    same = R["chain"][1:] == R["chain"][:-1]
    np.add.at(T, (st[:-1][same], st[1:][same]), 1.0)
    R["trans"] = T / T.sum(1, keepdims=True)
    np.savez_compressed(res_path(fold), **R)
    if verbose:
        print(f"fold {fold}: {ci} chains ({len(peps)} pep, {len(frs)} frag), {len(st)} residues", flush=True)
    # pieces
    allca = np.concatenate([np.asarray(q.ca, np.float32) for q in peps + frs])
    for L in LS:
        if os.path.exists(piece_path(fold, L)):
            continue
        start = np.where(R["pos"] + L <= R["clen"])[0].astype(np.int32)
        idx = start[:, None] + np.arange(L)[None]
        ca = np.zeros((len(start), L, 3), np.float32)
        for a in range(0, len(start), 20000):
            sl = idx[a:a + 20000]
            ca[a:a + 20000] = I.build_ca(R["phi"][sl].astype(float), R["psi"][sl].astype(float))
        np.savez_compressed(piece_path(fold, L), start=start, ca=ca, rca=allca[idx])
        if verbose:
            print(f"  L={L}: {len(start)} pieces", flush=True)


_RES, _PIECE = {}, {}


def residues(fold):
    if fold not in _RES:
        z = np.load(res_path(fold))
        _RES[fold] = {k: z[k] for k in z.files}
    return _RES[fold]


def pieces(fold, L):
    key = (fold, L)
    if key not in _PIECE:
        z = np.load(piece_path(fold, L))
        _PIECE[key] = {k: z[k] for k in z.files}
    return _PIECE[key]


def piece_torsions(fold, L, sel):
    """(phi, psi, codes) each (len(sel), L) for piece indices sel."""
    R = residues(fold); P = pieces(fold, L)
    idx = P["start"][np.asarray(sel, int)][:, None] + np.arange(L)[None]
    return R["phi"][idx].astype(float), R["psi"][idx].astype(float), R["code"][idx]


if __name__ == "__main__":
    print("free GB", I.free_gb())
    for f in range(5):
        build_fold(f)
