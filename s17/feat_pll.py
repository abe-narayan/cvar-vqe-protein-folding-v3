"""s17/feat_pll.py -- ESM-2 per-position amino-acid log-probabilities for the 126 targets.

WHY THIS AND NOT AN EMBEDDING.  A candidate window is a real fragment carrying its OWN source
sequence (`u["S"]`).  The question "does this fragment's sequence belong at this position of the
target" is a learned-sequence question that BLOSUM -- a position-INDEPENDENT substitution matrix,
already measured null in-band (SELECT s5.1, rho -0.013) -- cannot express.  ESM-2 answers it
directly: run the model on the TARGET sequence and read the per-position distribution over the
20 amino acids.  Scoring a candidate by

    -mean_r log P_esm( S_candidate[r] | target sequence )

is then a per-candidate feature built from a learned sequence representation, and it costs
**one forward pass per target** -- 126 in total -- rather than one per candidate window.

TWO CONDITIONING MODES, both cached, because they are different objects:
    `nat`     the model sees the whole target sequence and we read position r's logits.  Cheap
              (1 pass) but position r's distribution is contaminated by residue r itself.
    `mask`    residue r is replaced by <mask> before the pass -- the proper pseudo-likelihood.
              Costs n passes per target (~1,600 total), still small.

Both are stored.  The masked form is the honest one and is the arm reported; the unmasked form
is kept because a difference between them is itself diagnostic.

COMPUTE.  ESM-2-650M is ~2.6 GB of weights.  This module loads the model ONCE, runs, frees it,
and writes `s17/cache/feat_pll.npz`.  It is the only heavy step in this lane and is run alone.
Nothing here reads a structure, so it is fold-honest by construction (the model is a fixed
pretrained artefact, identical for every fold, exactly as the shipped distogram's ESM block is).
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

OUT = os.path.join(HERE, "cache", "feat_pll.npz")
MODEL_NAME = "esm2_t33_650M_UR50D"


def build(verbose=True):
    import torch
    import esm as esmlib
    torch.set_num_threads(3)

    tg = I.targets()
    model, alphabet = getattr(esmlib.pretrained, MODEL_NAME)()
    model.eval()
    bc = alphabet.get_batch_converter()
    aa_tok = np.array([alphabet.get_idx(a) for a in I.ALPHABET], int)
    mask_tok = alphabet.mask_idx

    seqs, LPn, LPm = [], [], []
    t0 = time.time()
    with torch.no_grad():
        for q, t in enumerate(tg):
            s = t["seq"]
            n = len(s)
            _, _, toks = bc([("x", s)])
            lg = model(toks)["logits"][0]                       # (n+2, V)
            lp = torch.log_softmax(lg, -1).numpy()[1:n + 1][:, aa_tok]
            LPn.append(lp.astype(np.float32))

            #: masked pseudo-likelihood: one pass per residue, batched in one go
            B = toks.repeat(n, 1)
            for r in range(n):
                B[r, r + 1] = mask_tok
            lgm = model(B)["logits"]
            lpm = torch.log_softmax(lgm, -1).numpy()
            lpm = np.stack([lpm[r, r + 1][aa_tok] for r in range(n)])
            LPm.append(lpm.astype(np.float32))
            seqs.append(s)
            if verbose and (q + 1) % 10 == 0:
                print(f"  pll {q+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
    np.savez_compressed(OUT, seqs=np.array(seqs, dtype=object),
                        lp_nat=np.array(LPn, dtype=object),
                        lp_mask=np.array(LPm, dtype=object))
    print(f"wrote {OUT}  ({time.time()-t0:.0f}s)", flush=True)


_P = None


def load():
    """seq -> (lp_nat (n,20), lp_mask (n,20)) log-probabilities over `I.ALPHABET`."""
    global _P
    if _P is None:
        if not os.path.exists(OUT):
            return {}
        z = np.load(OUT, allow_pickle=True)
        _P = {str(s): (np.asarray(a, np.float32), np.asarray(b, np.float32))
              for s, a, b in zip(z["seqs"], z["lp_nat"], z["lp_mask"])}
    return _P


if __name__ == "__main__":
    build()
