"""SHIFT agent -- the diagnostic that decides how to read `shift_model`'s 61.2 deg.

If the same model reaches TALOS-N-like accuracy on HELD-OUT FOLDED PROTEINS but 61 deg on
9-16-mers, the shortfall is DOMAIN SHIFT and the channel is intact.  If it reaches 60 deg on
held-out proteins too, the model is simply weak and the peptide number says nothing about the
channel.  These are opposite conclusions and the experiment costs minutes.

Also sweeps the context window (TALOS-N matches HEPTApeptides, +/-3; `shift_model` uses
+/-1) and the nucleus subsets our targets actually carry, so the cost of each missing
nucleus is measured rather than assumed (section 2.4 flagged it UNVERIFIED).

    python -m s14.shift_diag
"""
from __future__ import annotations

import os
import sys
import json
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import torch                                         # noqa: E402

from s13 import tors_common as T                     # noqa: E402
from s14.shift_bmrb import BACKBONE, RESULTS         # noqa: E402
from s14 import shift_model as SM                    # noqa: E402

torch.set_num_threads(2)
R2D = 180.0 / np.pi


def build_feats(rows, win, ref, glob, sc, keep_nuc=None):
    """Re-featurise the corpus with an arbitrary context window.  `rows` carry only
    +/-1 shift context, so wider windows use sequence context only beyond +/-1 -- which is
    itself the measurement: how much of TALOS-N's heptapeptide gain is SEQUENCE context."""
    nb = len(BACKBONE)
    nf_sh = 3 * nb * 2
    nf_aa = (2 * win + 1) * (len(SM.AA) + 1)
    X = np.zeros((len(rows), nf_sh + nf_aa), np.float32)
    for n, r in enumerate(rows):
        aas = [r["prev"], r["aa"], r["next"]]
        shs = [r["sh_prev"], r["sh"], r["sh_next"]]
        if keep_nuc is not None:
            shs = [{k: v for k, v in d.items() if k in keep_nuc} for d in shs]
        k = 0
        for w in (-1, 0, 1):
            j = 1 + w
            aa = aas[j]
            for a in BACKBONE:
                v = shs[j].get(a)
                if v is None:
                    X[n, k] = 0.0
                    X[n, k + 1] = 0.0
                else:
                    X[n, k] = float((v - ref.get((aa, a), glob[a])) / sc[a])
                    X[n, k + 1] = 1.0
                k += 2
        ctx = r.get("ctx")
        for off, w in enumerate(range(-win, win + 1)):
            if w == -1:
                aa = r["prev"]
            elif w == 0:
                aa = r["aa"]
            elif w == 1:
                aa = r["next"]
            elif ctx:
                aa = ctx.get(str(w), "-")
            else:
                aa = "-"
            X[n, k + SM.AAI.get(aa, len(SM.AA))] = 1.0
            k += len(SM.AA) + 1
    return X


def run(rows, ref, glob, sc, ycell, phi, psi, groups, keep_nuc, tag, ho):
    X = build_feats(rows, 1, ref, glob, sc, keep_nuc)
    tr = ~ho
    m = SM.fit(X[tr], ycell[tr], 324, "cls", groups=[g for g, k in zip(groups, tr) if k],
               tag=tag, verbose=False)
    P = SM.posterior(m, X[ho])
    d = T.decode(P)
    ep = T.wrap(d["phi"] - phi[ho]) * R2D
    es = T.wrap(d["psi"] - psi[ho]) * R2D
    sig = float(np.sqrt(((ep ** 2).mean() + (es ** 2).mean()) / 2))
    return {"tag": tag, "n_test": int(ho.sum()), "sigma": sig,
            "mae_phi": float(np.abs(ep).mean()), "mae_psi": float(np.abs(es).mean()),
            "gross_rate": float((np.hypot(ep, es) > 60).mean()),
            "cell_acc": float(np.mean(np.argmax(P, 1) == ycell[ho]))}


def main():
    t0 = time.time()
    rows = SM.load_corpus()
    ref, glob, sc = SM.reference(rows)
    _, cphi, cpsi = SM.corpus_xy(rows, ref, glob, sc)
    ycell = T.grid_bin(cphi, cpsi).astype(np.int64)
    groups = np.array([r["pdb"] for r in rows])
    prots = np.unique(groups)
    rng = np.random.default_rng(0)
    test = set(rng.permutation(prots)[:80].tolist())
    ho = np.isin(groups, list(test))
    print("corpus {} residues / {} proteins; held-out {} proteins, {} residues"
          .format(len(rows), len(prots), len(test), int(ho.sum())))
    print()
    print("{:34s} {:>7s} {:>8s} {:>8s} {:>8s} {:>8s}".format(
        "nucleus set", "sigma", "MAEphi", "MAEpsi", "gross", "cellacc"))
    out = {"n_residues": len(rows), "n_proteins": int(len(prots)),
           "n_heldout_proteins": len(test), "rows": []}
    sets = [
        ("ALL SIX  H+N+HA+CA+CB+C", tuple(BACKBONE)),
        ("H+N+HA+CA+CB  (no C-prime)", ("H", "N", "HA", "CA", "CB")),
        ("H+HA+CA+CB    (no N, no C)", ("H", "HA", "CA", "CB")),
        ("N+CA+CB+C     (no 1H)", ("N", "CA", "CB", "C")),
        ("CA+CB only", ("CA", "CB")),
        ("CA only", ("CA",)),
        ("NONE (sequence only)", ()),
    ]
    for lbl, ks in sets:
        r = run(rows, ref, glob, sc, ycell, cphi, cpsi, groups, set(ks), lbl, ho)
        out["rows"].append(r)
        print("{:34s} {:7.1f} {:8.1f} {:8.1f} {:8.3f} {:8.3f}".format(
            lbl, r["sigma"], r["mae_phi"], r["mae_psi"], r["gross_rate"], r["cell_acc"]))
    out["what"] = ("Held-out-PROTEIN accuracy of the same shift->torsion model that "
                   "s14/shift_model.py applies to the 9-16-mers. The gap between these "
                   "numbers and the peptide number is DOMAIN SHIFT.")
    out["secs"] = round(time.time() - t0, 1)
    with open(os.path.join(RESULTS, "shift_diag.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote s14/results/shift_diag.json  [{:.0f}s]".format(time.time() - t0))


if __name__ == "__main__":
    main()
