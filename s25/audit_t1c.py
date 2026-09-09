"""s25/audit_t1c.py -- THE DECISIVE EXPERIMENT L2 DID NOT RUN: location and width, separated.

s25/temper.py's SD(f) operator is NOT a width operator.  Measured (s25/audit_t1.py, B1) it moves
the posterior MEAN by -0.375 A at f=1.5 and -1.155 A at f=3.0, and the shift is separation-graded
(r = 0.977 with L1's own signed error across the five bands).  A Gaussian kernel in DISTANCE space,
row-normalised over a NON-UNIFORM bin grid (0.5 A spacing at 4 A, 3.5 A at 21 A), necessarily drags
mass toward the dense short-distance bins.  So the arm labelled "widen" was "widen AND contract by
~1 A", and L2's mechanism claim -- location matters, width is inert -- cannot be read off it.

THIS FILE RUNS THE TWO ARMS SEPARATELY.

    MPW(f)      MEAN-PRESERVING widening.  SD-convolve, then translate the mass back so the
                posterior's mean is exactly restored.  Width moves, location does not.
    MPWFIXW(f)  identical posterior, per-pair weight w frozen at the shipped value, so the
                width->weight channel (w = shell/(sd+0.5)^g, shell == 1, g == 1) is closed too.
    SHIFT(d)    pure translation of the shipped posterior by d Angstrom.  Location moves,
                width does not.

OPERATOR FORKS (rule 0; enumerated by the AUDIT lane, which has no stake in the outcome).

    functional     DECLARED the shipped Bayes risk via a genuine core.predict.Distogram, gated on
                   bit-exact reconstruction of the shipped risk table per target.  NOT TAKEN a
                   re-implementation, and NOT TAKEN altering the loss.
    basis          DECLARED point cloud, medoid frame, shipped uniform top-75 average -- identical
                   to temper.py so the two are directly comparable.  NOT TAKEN a built chain.
    readout        DECLARED per-target Ca-RMSD, paired against the shipped f=1 arm on the SAME
                   targets.  NOT TAKEN a pooled or re-ranked readout.
    normalisation  DECLARED translation by linear redistribution onto the shipped bin centres,
                   which preserves the first moment EXACTLY.  NOT TAKEN a nearest-bin roll, which
                   quantises the shift to the (non-uniform) bin spacing.
    null           DECLARED the shipped posterior, d = 0 / f = 1, an exact interior point of both
                   grids.  NOT TAKEN a best-of-grid oracle; the grids are reported in full so no
                   minimum is taken over them.
    THE LABEL      DECLARED Ca-RMSD AND, beside it, the ACHIEVED sd ratio and ACHIEVED mean shift
                   of every arm, so the parameter's name is checked against what it does -- the
                   exact audit that broke SD(f).  NOT TAKEN reporting the nominal parameter alone.

  Prediction registered BEFORE the run.  If L2's mechanism claim is right, MPW must be flat in f
  (width inert) and SHIFT must be steep in d (location live).  If instead MPW is steep, width is
  live and L2's mechanism is wrong.  If BOTH are flat past their MDEs, the honest reading is that
  the endpoint does not respond to the posterior at this amplitude at all, and L2's +0.0538 was
  the noise it already measures itself to be (0.38x its own MDE).
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I           # noqa: E402
from core import predict as dgm           # noqa: E402
from s24 import stats_lib as ST           # noqa: E402

TOPM = 75
FGRID = [1.0, 1.3, 1.5, 2.0, 3.0]
DGRID = [-0.9, -0.6, -0.3, 0.0, +0.3]


def _widen_sd(prob, C, f):
    if f <= 1.0:
        return prob.copy()
    m = (prob * C[None, :]).sum(1)
    v = (prob * (C[None, :] - m[:, None]) ** 2).sum(1)
    s = np.sqrt(np.maximum(v, 1e-12)) * np.sqrt(f * f - 1.0)
    d2 = (C[None, :] - C[:, None]) ** 2
    out = np.empty_like(prob)
    for p in range(len(prob)):
        K = np.exp(-d2 / (2.0 * max(s[p], 1e-6) ** 2))
        K /= K.sum(1, keepdims=True)
        q = prob[p] @ K
        out[p] = q / q.sum()
    return out


def _translate(prob, C, delta):
    """Move every bin's mass to C + delta and redistribute linearly onto C.

    Linear (mass-splitting) redistribution preserves the first moment exactly for any mass that
    stays inside [C[0], C[-1]]; mass pushed past an end is clamped there, which is reported via
    the ACHIEVED mean shift rather than assumed away.  `delta` may be a scalar or per-pair.
    """
    d = np.broadcast_to(np.asarray(delta, float).reshape(-1, 1), prob.shape)
    tgt = np.clip(C[None, :] + d, C[0], C[-1])
    k = np.clip(np.searchsorted(C, tgt, side="right") - 1, 0, len(C) - 2)
    lo, hi = C[k], C[k + 1]
    frac = (tgt - lo) / np.maximum(hi - lo, 1e-12)
    out = np.zeros_like(prob)
    r = np.repeat(np.arange(prob.shape[0]), prob.shape[1])
    np.add.at(out, (r, k.ravel()), (prob * (1.0 - frac)).ravel())
    np.add.at(out, (r, (k + 1).ravel()), (prob * frac).ravel())
    return out / out.sum(1, keepdims=True)


def _mean(prob, C):
    return (prob * C[None, :]).sum(1)


def _sd(prob, C):
    m = _mean(prob, C)
    return np.sqrt(np.maximum((prob * (C[None, :] - m[:, None]) ** 2).sum(1), 1e-12))


def _mpw(prob, C, f):
    """Widen by f, then translate each pair back so its mean is exactly restored."""
    if f <= 1.0:
        return prob.copy()
    m0 = _mean(prob, C)
    q = _widen_sd(prob, C, f)
    q = _translate(q, C, m0 - _mean(q, C))
    r = _translate(q, C, m0 - _mean(q, C))       # one refinement for end-clamped mass
    return r


def _avg(members):
    P = I.pairwise_rmsd(members); b = I.medoid(P)
    return I.superpose_batch(members, members[b]).mean(0)


def run():
    tg = I.targets()
    rows = []
    print("targets %d   f %s   delta %s" % (len(tg), FGRID, DGRID), flush=True)
    for c_i, t in enumerate(tg):
        pdb = t["pdb"]
        u = I.load_univ(pdb)
        seq = u["seq"]; nat = np.asarray(u["nat_ca"], float)
        pi, pj = I.pair_index(int(u["n"]))
        dg = I.distogram(pdb, seq, u["fold"])
        P0 = np.asarray(dg["prob"], float); C = np.asarray(dg["centres"], float)
        Wpool = np.asarray(u["W"], float)[I.pool_idx(u)]
        D = I.pair_dists(Wpool, pi, pj)

        d0 = dgm.Distogram(seq, P0, pi, pj)
        err = float(np.abs(np.asarray(dg["risk"], float) - np.asarray(d0._risk, float)).max())
        assert err < 1e-4, (pdb, err)
        w_real = np.asarray(d0.w, float)
        m0 = _mean(P0, C); s0 = _sd(P0, C)

        def emit(pr, fixw=False):
            d = dgm.Distogram(seq, pr, pi, pj)
            rk = np.asarray(d._risk, float)
            if fixw:
                rk = rk / np.maximum(np.asarray(d.w, float)[:, None], 1e-12) * w_real[:, None]
            sc = np.asarray(I.shipped_score({"grid": d.grid, "risk": rk}, D), float)
            return float(I.ca_rmsd(_avg(Wpool[np.argsort(sc, kind="stable")[:TOPM]]), nat))

        row = {"pdb": pdb, "n": int(u["n"]), "fold": int(u["fold"]), "recon_err": err,
               "MPW": [], "MPWFIXW": [], "SHIFT": [],
               "mpw_sdratio": [], "mpw_dmean": [], "sh_sdratio": [], "sh_dmean": []}
        for f in FGRID:
            pr = _mpw(P0, C, f)
            row["MPW"].append(emit(pr))
            row["MPWFIXW"].append(emit(pr, fixw=True))
            row["mpw_sdratio"].append(float((_sd(pr, C) / s0).mean()))
            row["mpw_dmean"].append(float((_mean(pr, C) - m0).mean()))
        for dl in DGRID:
            pr = _translate(P0, C, dl) if dl != 0.0 else P0.copy()
            row["SHIFT"].append(emit(pr))
            row["sh_sdratio"].append(float((_sd(pr, C) / s0).mean()))
            row["sh_dmean"].append(float((_mean(pr, C) - m0).mean()))
        rows.append(row)
        if (c_i + 1) % 20 == 0:
            print("  %d/%d" % (c_i + 1, len(tg)), flush=True)
            ST.save_atomic(os.path.join(RES, "audit_t1c.json"),
                           {"rows": rows, "fgrid": FGRID, "dgrid": DGRID, "topm": TOPM},
                           rows=rows, n_expected=len(tg), module_file=__file__)
    ST.save_atomic(os.path.join(RES, "audit_t1c.json"),
                   {"rows": rows, "fgrid": FGRID, "dgrid": DGRID, "topm": TOPM},
                   complete_keys=("MPW", "MPWFIXW", "SHIFT", "recon_err"),
                   rows=rows, n_expected=len(tg), module_file=__file__)
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(os.path.join(RES, "audit_t1c.json")))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    pdbs = [r["pdb"] for r in rows]
    M = np.array([r["MPW"] for r in rows], float)
    MF = np.array([r["MPWFIXW"] for r in rows], float)
    S = np.array([r["SHIFT"] for r in rows], float)
    base = M[:, 0]
    b0 = DGRID.index(0.0)
    assert abs(S[:, b0].mean() - base.mean()) < 1e-9, "the two identity arms must agree exactly"

    print("\nn = %d.  Point cloud, shipped uniform top-75 average; only the SELECTION varies." % len(rows))
    print("  risk-table reconstruction: max err %.2e" % max(r["recon_err"] for r in rows))
    print("\n  MEAN-PRESERVING WIDENING -- width moves, location does not")
    print("  %8s%12s%12s%12s%12s%14s" % ("f", "RMSD", "RMSD fixw", "sd ratio", "d(mean)", "vs f=1"))
    for a, f in enumerate(FGRID):
        print("  %8.2f%12.4f%12.4f%12.3f%+12.4f%+14.4f"
              % (f, M[:, a].mean(), MF[:, a].mean(),
                 np.mean([r["mpw_sdratio"][a] for r in rows]),
                 np.mean([r["mpw_dmean"][a] for r in rows]), (M[:, a] - base).mean()))
    print("\n  PURE TRANSLATION -- location moves, width does not")
    print("  %8s%12s%12s%12s%14s" % ("delta", "RMSD", "sd ratio", "d(mean)", "vs delta=0"))
    for a, d in enumerate(DGRID):
        print("  %+8.2f%12.4f%12.3f%+12.4f%+14.4f"
              % (d, S[:, a].mean(), np.mean([r["sh_sdratio"][a] for r in rows]),
                 np.mean([r["sh_dmean"][a] for r in rows]), (S[:, a] - base).mean()))

    print("\n  THROUGH stats_lib (negative = better than the shipped posterior):")
    for a, f in enumerate(FGRID):
        if f == 1.0:
            continue
        print(ST.fmt(ST.compare(M[:, a], base, fold, names=pdbs,
                                label="MPW f=%.2f (pure width) vs shipped" % f)))
    for a, d in enumerate(DGRID):
        if d == 0.0:
            continue
        print(ST.fmt(ST.compare(S[:, a], base, fold, names=pdbs,
                                label="SHIFT d=%+.2f A (pure location) vs shipped" % d)))


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
