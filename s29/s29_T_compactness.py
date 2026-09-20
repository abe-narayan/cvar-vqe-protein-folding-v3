#!/usr/bin/env python
"""s29/s29_T_compactness.py -- LANE T: IS THE PHYSICS/REALISM FAMILY ON THE COMPACTNESS AXIS?

WHY THIS FILE EXISTS, AND WHY IT IS NOT THE MEASUREMENT THAT WAS ASKED FOR. The coordinator asked
for S8's free-energy stage (`strain`, `E_free`, `F_qh`, `F_boltz`, `width`, `S_msf`) to be
correlated against compactness, on the 1 of 24 targets it completed. **That stage does not exist.**
`s8/relax.py`, `s8/test_relax.py`, `s8/relax_fe.json`, `s8/relax_best.json`, `s8/relax_sweep.json`,
`s8/relax_report.txt` and `s8/relax_findings.md` are absent from disk AND from git history
(`git log --all -- <path>` returns nothing for all seven). Only the prose merged into
`docs/FINDINGS.md` section B survives, so "committed and resumable as `python -m s8.relax best`"
is not true and no `width` or `S_msf` column exists anywhere to correlate.

WHAT THIS MEASURES INSTEAD, on artefacts that do exist. Lane L's objection is that a
native-free physics score at peptide length is plausibly compactness-loaded, in which case the
free-energy class is not orthogonal to the realism axis. That question is answerable for EVERY
channel the project actually owns, on real pools, without OpenMM and without a native in the
selection half:

    for each of the 32 S27 channels, over the 500 members of each pool:
      rho_rg     = Spearman(channel, Rg of the member)                       NATIVE-FREE
      rho_inband = Spearman(channel, member's CA-RMSD to the native), rr < 3 A     ORACLE
      rho_part   = the same, PARTIALLING OUT Rg (Spearman residuals)              ORACLE

`rho_rg` says whether the channel IS compactness; `rho_inband - rho_part` says how much of its
in-band skill is compactness. The ORACLE columns are diagnostics and choose nothing.

    python s29/s29_T_compactness.py [--limit N]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

OUT = os.path.join(HERE, "results", "s29_T_compactness.json")


def rg_of(W):
    W = np.asarray(W, float)
    c = W.mean(1, keepdims=True)
    return np.sqrt(((W - c) ** 2).sum(-1).mean(1))


def _rank(x):
    from scipy.stats import rankdata
    return rankdata(np.asarray(x, float))


def partial_spearman(x, y, z):
    """Spearman correlation of x and y with z partialled out, on ranks (Pearson of the residuals
    of the rank vectors, the standard Spearman-partial)."""
    rx, ry, rz = _rank(x), _rank(y), _rank(z)
    def resid(a, b):
        b = b - b.mean()
        return a - a.mean() - (np.dot(a - a.mean(), b) / max(np.dot(b, b), 1e-30)) * b
    ex, ey = resid(rx, rz), resid(ry, rz)
    den = np.sqrt(np.dot(ex, ex) * np.dot(ey, ey))
    return float(np.dot(ex, ey) / den) if den > 0 else float("nan")


def target_row(pdb):
    from scipy.stats import spearmanr
    from s27 import run_pool as RP
    cand, ch, _ = RP.channels_for(pdb)
    rr = np.asarray(cand.oracle_rr, float)              # ORACLE, diagnostic only
    rg = rg_of(cand.W)                                   # native-free
    inb = rr < 3.0
    row = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), k=int(cand.k),
               n_inband=int(inb.sum()), rg_mean=float(rg.mean()),
               rho_rg_rr=float(spearmanr(rg, rr).correlation),
               rho_rg_rr_inband=float(spearmanr(rg[inb], rr[inb]).correlation) if inb.sum() > 5 else float("nan"),
               ch={})
    for name, v in ch.items():
        v = np.asarray(v, float)
        if not np.isfinite(v).all() or np.std(v) == 0:
            fin = np.isfinite(v)
            if fin.sum() < 50 or np.std(v[fin]) == 0:
                continue
            v = np.where(fin, v, np.nanmax(v[fin]) + 1.0)
        cell = dict(rho_rg=float(spearmanr(v, rg).correlation))
        if inb.sum() > 10 and np.std(v[inb]) > 0:
            cell["rho_inband"] = float(spearmanr(v[inb], rr[inb]).correlation)
            cell["rho_inband_partial_rg"] = partial_spearman(v[inb], rr[inb], rg[inb])
        row["ch"][name] = cell
    return row


def summarise(rows):
    from s24 import stats_lib as ST
    names = sorted({k for r in rows for k in r["ch"]})
    folds = ST.pinned_folds([r["pdb"] for r in rows])
    out = dict(kind="s29_T_compactness", lane="T", n=len(rows),
               note=("the S8 free-energy stage is ABSENT from disk and git; this measures the "
                     "compactness loading of the channels that do exist"),
               rho_rg_rr_median=float(np.median([r["rho_rg_rr"] for r in rows])),
               rho_rg_rr_inband_median=float(np.nanmedian([r["rho_rg_rr_inband"] for r in rows])),
               channels={})
    for nm in names:
        g = np.array([r["ch"][nm]["rho_rg"] for r in rows if nm in r["ch"]], float)
        ib = np.array([r["ch"][nm].get("rho_inband", np.nan) for r in rows if nm in r["ch"]], float)
        pa = np.array([r["ch"][nm].get("rho_inband_partial_rg", np.nan) for r in rows if nm in r["ch"]], float)
        fo = folds[[i for i, r in enumerate(rows) if nm in r["ch"]]]
        ok = np.isfinite(g)
        cell = dict(n_targets=int(ok.sum()),
                    rho_rg_median=float(np.median(g[ok])),
                    rho_rg_abs_median=float(np.median(np.abs(g[ok]))),
                    fold_ci_rho_rg=ST.compare(g[ok], np.zeros(int(ok.sum())), fo[ok],
                                              label="rho_rg " + nm, seed_parts=("s29Tcp",))["ci95_fold"])
        if np.isfinite(ib).sum() > 30:
            m = np.isfinite(ib) & np.isfinite(pa)
            cell.update(rho_inband_median=float(np.median(ib[m])),
                        rho_inband_partial_median=float(np.median(pa[m])),
                        share_of_inband_from_rg=float(1.0 - np.median(np.abs(pa[m])) /
                                                      max(np.median(np.abs(ib[m])), 1e-9)),
                        n_inband_targets=int(m.sum()))
        out["channels"][nm] = cell
    return out


def main(argv=None):
    from s24 import stats_lib as ST
    from s12 import instrument as I
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(argv)
    pdbs = [t["pdb"] for t in I.targets()]
    if a.limit:
        pdbs = pdbs[:a.limit]
    rows = []
    for q, pdb in enumerate(pdbs):
        rows.append(target_row(pdb))
        if (q + 1) % 30 == 0:
            print(f"  [{q+1}/{len(pdbs)}]", flush=True)
    out = summarise(rows)
    out["rows"] = rows
    ST.save_atomic(OUT, out)
    print(f"n = {out['n']}; Rg's own in-band skill (ORACLE): rho {out['rho_rg_rr_inband_median']:+.3f}")
    print(f"{'channel':22s} {'rho(X,Rg)':>10s} {'|rho|':>7s} {'in-band':>9s} {'partial':>9s} {'Rg share':>9s}")
    for nm, c in sorted(out["channels"].items(), key=lambda kv: -abs(kv[1]["rho_rg_median"])):
        print(f"{nm:22s} {c['rho_rg_median']:+10.3f} {c['rho_rg_abs_median']:7.3f} "
              f"{c.get('rho_inband_median', float('nan')):+9.3f} "
              f"{c.get('rho_inband_partial_median', float('nan')):+9.3f} "
              f"{c.get('share_of_inband_from_rg', float('nan')):9.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
