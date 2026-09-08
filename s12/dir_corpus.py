"""D3a. Build the SIGN-LABEL training corpus: a leave-fold-out distogram for every one of
the 787 peptides in `peptide_db`, not just the 126 tuning targets.

Discipline: peptide P (fold f) is scored by `fold_model(f)`, which never saw fold f.  That
is the same leave-fold-out rule the shipped pipeline uses for the tuning targets, so the
error signs recorded here are the errors a deployed distogram really makes.  Downstream,
the sign head that scores target T is trained WITHOUT fold(T), so T's own errors never
train the head that judges T.

ESM: served from `s12/cache/dir_esm_pep.npz` (787 peptides, extracted from the compact
`s12/esm_bank.py` bank).  `esm_cache.npz` is NEVER loaded.  `core.data.esm_embed` returns
the same whitened PCA-32 the pipeline uses (esm_pca.npz), so this is the shipped feature
path; the only difference is float16 storage in the bank, validated below against the
126 cached shipped distograms.

Writes s12/cache/dir_corpus.npz.
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
import torch
torch.set_num_threads(2)

from s12 import instrument as I

OUT = os.path.join(ROOT, "s12", "cache", "dir_corpus.npz")
PEPESM = os.path.join(ROOT, "s12", "cache", "dir_esm_pep.npz")


def install_esm_bank():
    """Serve esm_embed / esm_contacts for the peptide corpus from the compact bank."""
    from core import data as cdata
    z = np.load(PEPESM, allow_pickle=True)
    bank = {str(s): (np.asarray(a, np.float32), np.asarray(c, np.float32))
            for s, a, c in zip(z["seqs"], z["pca32"], z["con"])}
    _emb, _con = cdata.esm_embed, cdata.esm_contacts
    hot = cdata._hot()

    def esm_embed(sequence):
        if sequence in hot:
            return _emb(sequence)
        if sequence in bank:
            return bank[sequence][0]
        raise RuntimeError(f"no ESM for {sequence!r}")

    def esm_contacts(sequence):
        if sequence in hot:
            return _con(sequence)
        if sequence in bank:
            return bank[sequence][1]
        raise RuntimeError(f"no ESM for {sequence!r}")

    cdata.esm_embed = esm_embed
    cdata.esm_contacts = esm_contacts
    return bank


def main():
    from core import pipeline as pl
    from core import predict as dgm
    recs = np.load(os.path.join(ROOT, "peptide_db.npz"), allow_pickle=True)["records"]
    folds = json.load(open(os.path.join(ROOT, "peptide_folds.json")))["assign"]
    install_esm_bank()
    pl.guard_esm()

    tg = {t["pdb"]: t for t in I.targets()}
    out = {}
    keys, t0 = [], time.time()
    val = []
    for k, r in enumerate(recs):
        seq, pdb = r["seq"], r["pdb"]
        n = len(seq)
        if n < 8 or seq not in folds:
            continue
        fold = int(folds[seq])
        ca = np.asarray(r["ca"], float)
        if ca.shape[0] != n or not np.isfinite(ca).all():
            continue
        i, j = I.pair_index(n)
        if i.size == 0:
            continue
        d = dgm.Distogram.for_target(seq, model=pl.fold_model(fold))
        assert np.array_equal(d.i, i) and np.array_equal(d.j, j), pdb
        dtrue = np.linalg.norm(ca[i] - ca[j], axis=-1)
        out[f"{pdb}/prob"] = np.asarray(d.prob, np.float16)
        out[f"{pdb}/exp"] = np.asarray(d.expected, np.float32)
        out[f"{pdb}/sd"] = np.asarray(d.sd, np.float32)
        out[f"{pdb}/dtrue"] = dtrue.astype(np.float32)
        keys.append((pdb, seq, n, fold))
        if pdb in tg:                       # validation vs the cached shipped distogram
            dg = I.distogram(pdb)
            val.append(float(np.abs(np.asarray(dg["expected"], float) - d.expected).max()))
        if (k + 1) % 50 == 0:
            print(f"  {k+1}/{len(recs)} n_ok={len(keys)} [{time.time()-t0:.0f}s "
                  f"free={I.free_gb():.1f}]", flush=True)
    out["keys"] = np.array([f"{a}\t{b}\t{c}\t{dd}" for a, b, c, dd in keys], dtype=object)
    np.savez_compressed(OUT, **out)
    print("wrote", OUT, len(keys), "peptides", flush=True)
    print(f"validation vs shipped distogram on {len(val)} tuning targets: "
          f"max|dexpected| mean={np.mean(val):.4f} max={np.max(val):.4f}", flush=True)
    I.write("dir_corpus_meta", dict(n_peptides=len(keys),
                                    n_pairs=int(sum(len(out[f'{a}/exp']) for a, *_ in keys)),
                                    val_max_abs_dexp_mean=float(np.mean(val)),
                                    val_max_abs_dexp_max=float(np.max(val)),
                                    n_validated=len(val)))


if __name__ == "__main__":
    main()
