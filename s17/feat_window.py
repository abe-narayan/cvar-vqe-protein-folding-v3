"""s17/feat_window.py -- F6: PER-CANDIDATE ESM EMBEDDINGS, computed for the BAND only.

THE COST GATE, declared in the pre-registration before any of this ran.  Per-candidate ESM
embeddings DO NOT EXIST in this repo:

    esm_cache.npz is 1.5 GB keyed by FULL library sequence; the universes hold length-n WINDOWS
    (`u["S"]`, int8 codes) with no window->parent map.  Featurising the full universe is
    13k-27k windows x 126 targets ~ 1.6M short sequences -- out of budget on a box at 100% CPU
    with ~11 GB of RAM shared between four agents.

But the question is an IN-BAND question, so only the band needs embedding:

    shipped top-25 of K = 500   ->  <= 3,150 sequences  (2,300-3,000 unique)
    shipped top-75 of K = 500   ->  <= 9,450 sequences  (6,000-8,000 unique)

which is affordable.  This module featurises exactly those and nothing else, so no compute is
spent on candidates whose ranking was never in question.

WHAT IT MEASURES.  Each candidate window is a real fragment carrying its own source sequence of
the same length n as the target.  So the target's per-residue embedding matrix and the window's
are ALIGNED BY CONSTRUCTION, and "how similar are these two sequences in the language model's
space" is a per-residue quantity, not a pooled one.

    emb_cos     -mean_r cos( E_target[r], E_window[r] )      residue-aligned agreement
    emb_l2       mean_r || Z_target[r] - Z_window[r] ||       in the shipped whitened PCA-32 space
    emb_pool     -cos( mean_r E_target[r], mean_r E_window[r] )   the pooled (weaker) form

CONTROLS.  `blosum_sim` -- the position-independent substitution matrix, i.e. the classical
answer to the same question, already measured null in-band -- plus the constant alpha-helix and
beta-strand, and random-in-band.  If `emb_*` does not separate from `blosum_sim` it is BLOSUM
re-derived in 1280 dimensions, not a representation result.

MODEL CHOICE is a declared parameter, not a silent one.  `esm2_t33_650M_UR50D` is the shipped
model; `esm2_t6_8M_UR50D` is a 30 MB scout that costs ~80x less and is used ONLY when the box
cannot hold the large model.  Whichever ran is recorded in the artefact and quoted with every
number, because they are different instruments.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "3")
os.environ.setdefault("MKL_NUM_THREADS", "3")

from s12 import instrument as I            # noqa: E402
from s17 import sel_lib as L               # noqa: E402

BIG = "esm2_t33_650M_UR50D"
SMALL = "esm2_t6_8M_UR50D"


def out_path(model, B):
    return os.path.join(HERE, "cache", f"feat_win_{model}_B{B}.npz")


def band_sequences(B=75, K=500, verbose=True, targets=None):
    """The unique source sequences of the shipped objective's own top-B, per target."""
    tg = targets if targets is not None else I.targets()
    per, uniq = {}, {}
    for t in tg:
        p = L.pack(t["pdb"], K, want=("D",))
        sc = L.shipped(p)
        rk = np.argsort(sc, kind="stable")[:B]
        S = np.asarray(I.load_univ(t["pdb"])["S"], np.int8)[p["idx"]][rk]
        seqs = [I.decode(r) for r in S]
        per[t["pdb"]] = (rk.tolist(), seqs)
        for s in seqs:
            uniq[s] = None
    if verbose:
        print(f"  band B={B}: {sum(len(v[1]) for v in per.values())} windows, "
              f"{len(uniq)} unique sequences", flush=True)
    return per, list(uniq.keys())


def embed(seqs, model=BIG, batch=16, verbose=True):
    """seq -> (n, D) final-layer per-residue representation."""
    import torch
    import esm as esmlib
    torch.set_num_threads(3)
    mdl, alphabet = getattr(esmlib.pretrained, model)()
    mdl.eval()
    bc = alphabet.get_batch_converter()
    layer = mdl.num_layers
    out = {}
    t0 = time.time()
    by_len = {}
    for s in seqs:
        by_len.setdefault(len(s), []).append(s)
    done = 0
    with torch.no_grad():
        for n, group in sorted(by_len.items()):
            for a in range(0, len(group), batch):
                ch = group[a:a + batch]
                _, _, toks = bc([(f"w{q}", s) for q, s in enumerate(ch)])
                rep = mdl(toks, repr_layers=[layer])["representations"][layer].numpy()
                for q, s in enumerate(ch):
                    out[s] = rep[q, 1:n + 1].astype(np.float16)
                done += len(ch)
                if verbose and done % 400 < batch:
                    print(f"  embed {done}/{len(seqs)}  ({time.time()-t0:.0f}s)", flush=True)
    return out


def build(B=75, model=BIG, verbose=True, targets=None):
    per, seqs = band_sequences(B, verbose=verbose, targets=targets)
    emb = embed(seqs, model=model, verbose=verbose)
    keys = list(emb.keys())
    np.savez_compressed(out_path(model, B), model=model, B=B,
                        seqs=np.array(keys, dtype=object),
                        emb=np.array([emb[s] for s in keys], dtype=object),
                        pdbs=np.array(list(per.keys()), dtype=object),
                        rk=np.array([per[k][0] for k in per], dtype=object),
                        wseq=np.array([per[k][1] for k in per], dtype=object))
    print(f"wrote {out_path(model, B)}", flush=True)


def load(B=75, model=BIG):
    p = out_path(model, B)
    if not os.path.exists(p):
        return None
    z = np.load(p, allow_pickle=True)
    emb = {str(s): np.asarray(v, np.float32) for s, v in zip(z["seqs"], z["emb"])}
    band = {str(k): (np.asarray(r, int), [str(x) for x in w])
            for k, r, w in zip(z["pdbs"], z["rk"], z["wseq"])}
    return dict(model=str(z["model"]), B=int(z["B"]), emb=emb, band=band)


if __name__ == "__main__":
    b = int(sys.argv[1]) if len(sys.argv) > 1 else 75
    m = sys.argv[2] if len(sys.argv) > 2 else BIG
    build(b, m)
