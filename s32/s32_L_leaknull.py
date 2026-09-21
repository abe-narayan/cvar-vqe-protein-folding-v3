"""LANE L -- the leakage filter must be tested against ITS OWN NULL before it is used.

The first `long40` admission pass rejected 140 of 170 monomers for reaching identity >= 0.4
to a member of the peptide/fragment banks under the window-level (`norm="shorter"`)
normalisation, leaving ZERO targets.  Project memory
(`containment-threshold-is-at-the-null`) records exactly this failure at 9-16 residues:
*random sequences score 0.56-0.63 against the fragment bank*, so a containment threshold
measures chance, not leakage.

A 40-60 residue target scored against a bank whose members are 9-25 residues, normalised by
the SHORTER sequence, is the same construction at worse odds: 4 matches out of a 9-residue
fragment is already 0.44.

So this module measures the null the filter is supposed to beat.  For each candidate target
it computes the filter's statistic on (a) the real sequence and (b) shuffled sequences with
the SAME amino-acid composition and the same length -- a zero-information control that is
PLAUSIBLE rather than uniform (memory: `zero-information-control-must-be-plausible`).  A
threshold at which the real rejection rate is not far above the shuffled rate is at the
null and is not a leakage filter.

Three candidate statistics are compared at once, because the choice between them IS the
result:

    shorter   NW identity normalised by the shorter sequence  (what failed)
    longer    NW identity normalised by the longer sequence   (the pinned convention)
    verbatim  the longest bank member appearing as an exact substring of the target,
              as a fraction of the target length

    python -m s32.s32_L_leaknull run
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "s32", "results")

N_SHUF = 3
THRESH = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)


def run(n_targets=60, seed=0, verbose=True):
    from core import data as cdata
    import peptide_db as pdb
    import fragment_db as fdb

    rng = np.random.default_rng(seed)
    cen = json.load(open(os.path.join(RESULTS, "L1_census.json")))
    cand = sorted([r for r in cen["rows"] if 40 <= r["n"] <= 60],
                  key=lambda r: r["pdb"])[:n_targets]
    bank = [p.seq for p in pdb.load()] + [f.seq for f in fdb.load()]
    bank_set = set(bank)
    if verbose:
        print(f"{len(cand)} targets vs {len(bank)} bank members", flush=True)

    def stats(seq):
        sh = float(cdata.identity_many(seq, bank, norm="shorter").max())
        lo = float(cdata.identity_many(seq, bank).max())
        vb = max((len(b) for b in bank_set if b in seq), default=0) / len(seq)
        return sh, lo, vb

    real, shuf = [], []
    for k, r in enumerate(cand):
        real.append(stats(r["seq"]))
        chars = list(r["seq"])
        for _ in range(N_SHUF):
            rng.shuffle(chars)
            shuf.append(stats("".join(chars)))
        if verbose and (k + 1) % 10 == 0:
            print(f"  {k+1}/{len(cand)}", flush=True)

    real = np.array(real); shuf = np.array(shuf)
    names = ["shorter", "longer", "verbatim"]
    out = {"n_targets": len(cand), "n_shuffles_per_target": N_SHUF,
           "n_bank": len(bank), "thresholds": list(THRESH), "stats": {}}
    for a, nm in enumerate(names):
        rows = {"real_mean": float(real[:, a].mean()),
                "real_max": float(real[:, a].max()),
                "shuffled_mean": float(shuf[:, a].mean()),
                "shuffled_max": float(shuf[:, a].max()),
                "reject_rate": {}}
        for t in THRESH:
            rows["reject_rate"][str(t)] = {
                "real": float((real[:, a] >= t).mean()),
                "shuffled": float((shuf[:, a] >= t).mean())}
        out["stats"][nm] = rows
    path = os.path.join(RESULTS, "L1c_leak_null.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        for nm in names:
            s = out["stats"][nm]
            print(f"\n{nm}:  real mean {s['real_mean']:.3f} (max {s['real_max']:.3f})   "
                  f"SHUFFLED mean {s['shuffled_mean']:.3f} (max {s['shuffled_max']:.3f})")
            print("   threshold   reject(real)  reject(shuffled)")
            for t in THRESH:
                d = s["reject_rate"][str(t)]
                print(f"   {t:9.2f}   {d['real']:11.3f}   {d['shuffled']:14.3f}")
        print("\nwrote", path, flush=True)
    return out


if __name__ == "__main__":
    run()
