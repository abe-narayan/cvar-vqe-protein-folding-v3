#!/usr/bin/env python
"""s27/run_wave3.py -- PREREG addendum 2 (H13): SS_MATCH2 (H-bond secondary-structure call vs
universe-fitted propensities), the SS_MATCH weight grid completed at 0.75, and the fixed
propensity-contrast gate.  Point cloud, top-75, rows in `s27/results/wave3_rows.jsonl`.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import geometry as geo           # noqa: E402
from s12 import instrument as I            # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import ham_lib as HL              # noqa: E402
from s27 import run_pool as RP             # noqa: E402

ROWS = os.path.join(RP.RESULTS, "wave3_rows.jsonl")


def hbond_matrix(bb, thresh=-0.5):
    """(k, n, n) boolean: donor N-H of residue a bonded to acceptor C=O of residue b."""
    N, C, O = bb["N"], bb["C"], bb["O"]
    k, n = N.shape[:2]
    u = C[:, :-1] - O[:, :-1]
    u = u / np.maximum(np.linalg.norm(u, axis=-1, keepdims=True), 1e-9)
    Hh = np.full_like(N, np.nan); Hh[:, 1:] = N[:, 1:] + u
    q = 0.084 * 332.0
    B = np.zeros((k, n, n), bool)
    for a in range(1, n):
        for b in range(n):
            if abs(a - b) < 2:
                continue
            rON = np.linalg.norm(O[:, b] - N[:, a], axis=-1); rCH = np.linalg.norm(C[:, b] - Hh[:, a], axis=-1)
            rOH = np.linalg.norm(O[:, b] - Hh[:, a], axis=-1); rCN = np.linalg.norm(C[:, b] - N[:, a], axis=-1)
            e = q * (1 / np.maximum(rON, .5) + 1 / np.maximum(rCH, .5) - 1 / np.maximum(rOH, .5) - 1 / np.maximum(rCN, .5))
            B[:, a, b] = e < thresh
    return B


def ss_call(B):
    """Per-residue helix / strand flags from the H-bond matrix (DSSP-like minimal rules)."""
    k, n = B.shape[:2]
    helix = np.zeros((k, n), bool); strand = np.zeros((k, n), bool)
    for i in range(n - 5):
        turn = B[:, i + 4, i] & B[:, i + 5, i + 1]          # two consecutive i -> i+4 bonds
        for s in range(1, 5):
            helix[:, i + s] |= turn
    for a in range(n):
        for b in range(n):
            if abs(a - b) >= 3:
                strand[:, a] |= B[:, a, b]; strand[:, b] |= B[:, a, b]
    strand &= ~helix
    return helix, strand


def main():
    done = set()
    if os.path.exists(ROWS):
        for line in open(ROWS, encoding="utf-8"):
            try:
                done.add(json.loads(line)["pdb"])
            except Exception:
                pass
    t0 = time.time()
    pdbs = P.targets()
    for nn, pdb in enumerate(pdbs):
        if pdb in done:
            continue
        cand, ch, _ = RP.channels_for(pdb)
        k = cand.k; key = RP.rng_for(pdb, "tiekey").random(k)
        u = I.load_univ(pdb)
        # universe-fitted per-residue-type helix / strand propensity from the H-bond SS call
        Ub = geo.build_backbone_batch(np.asarray(u["PHI"], float)[:3000], np.asarray(u["PSI"], float)[:3000])
        Ub = {kk: np.asarray(v, float) for kk, v in Ub.items()}
        hU, sU = ss_call(hbond_matrix(Ub))
        S = np.asarray(u["S"], int)[:3000]
        ph = np.zeros(20); ps = np.zeros(20); cnt = np.zeros(20)
        np.add.at(ph, S.ravel(), hU.ravel()); np.add.at(ps, S.ravel(), sU.ravel()); np.add.at(cnt, S.ravel(), 1.0)
        ph = (ph + 1) / (cnt + 2); ps = (ps + 1) / (cnt + 2)
        aa = HL.codes(cand.seq)
        bb = {kk: np.asarray(v, float) for kk, v in geo.build_backbone_batch(cand.PHI, cand.PSI).items()}
        hC, sC = ss_call(hbond_matrix(bb))
        # energy: mismatch = helix where the residue type is rarely helical, strand where rarely strand
        E2 = (hC * (1 - ph[aa])[None, :] + sC * (1 - ps[aa])[None, :]).sum(1)
        ch2 = dict(ch); ch2["SS_MATCH2"] = E2
        rows = []
        base = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold))

        def rm(idx):
            return RP.oracle_rmsd_of_set(cand, idx)

        rows.append(dict(base, config="DIS", rmsd=rm(RP.topm(ch["DIS"], 75, key))))
        rows.append(dict(base, config="SS_MATCH2", rmsd=rm(RP.topm(E2, 75, key))))
        for w in (0.25, 0.5, 0.75, 1.0):
            rows.append(dict(base, config=f"DIS+{w}*SS_MATCH2", rmsd=rm(RP.topm(RP.combine(ch2, ["DIS", "SS_MATCH2"], [1, w]), 75, key))))
        rows.append(dict(base, config="DIS+0.75*SS_MATCH", rmsd=rm(RP.topm(RP.combine(ch, ["DIS", "SS_MATCH"], [1, 0.75]), 75, key))))
        # fixed contrast gate
        c_t = float(np.mean([abs(HL.CF_A[a] - HL.CF_B[a]) for a in cand.seq]))
        rows.append(dict(base, config="DIS+gate(SS_MATCH)", w_t=c_t, rmsd=rm(RP.topm(RP.combine(ch, ["DIS", "SS_MATCH"], [1, 0.5 * c_t / 0.66]), 75, key))))
        rows.append(dict(base, config="DIS+SS_MATCH+SS_MATCH2", rmsd=rm(RP.topm(RP.combine(ch2, ["DIS", "SS_MATCH", "SS_MATCH2"], [1, 0.5, 0.5]), 75, key))))
        with open(ROWS, "a", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        print(f"  [{nn+1}/{len(pdbs)}] {pdb} ({(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", ROWS)


if __name__ == "__main__":
    main()
