"""SPRINT 14, coordinator -- does the STRUCTURAL objective break the discrimination floor?

THE BOTTLENECK, and why this is the last question worth asking.

The ENER workstream established that discrimination is a function of the QUALITY GAP between
two structures rather than their separation, and found a hard floor:

    no PHYSICAL objective exceeds 0.511 pairwise accuracy when two structures differ by
    less than 0.25 A.  leg_steric needs a 0.58 A gap for 55%, Legacy 1.31 A, AMBER 1.58 A.

That single threshold explains everything else this sprint measured: the ~1.9 A selection gap
that reproduces on two disjoint target sets, its survival at the certified optimum, the
collapse of coordinate descent from 1.295 A to 3.68 A, and the budget saturation.

But their sweep covered the physical objectives and the Legacy components. **It never
included the native-free STRUCTURAL objective this sprint built** -- the one whose certified
global optimum is 0.885 A better than random where Legacy's is 0.139 A worse. That objective
is better on every axis measured so far, so whether it also breaks the floor is the single
most decision-relevant question left: if it does, the bottleneck is not universal and there
is a route; if it does not, the floor is a property of every objective available and the
sprint's negative is complete.

This runs on the FULL 262,144-configuration enumeration of all nine enumerated targets --
no sampling, no proposal, so none of the proposal confound that forced the C7 retraction.

Uses the ENER workstream's own `pair_accuracy` so the numbers compose with theirs.

ORACLE for scoring only; the objective itself is native-free.

Run:
    python -m s14.discrim
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402
from s14 import ener_lib as E              # noqa: E402
from s14 import hamil as H                 # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(RESULTS, exist_ok=True)

ENUM = ["1CS9", "2MK7", "2P5H", "6EY3", "6F3V", "6S0N", "7N2I", "8IS3", "9UV5"]
GAPS = [0.0, 0.10, 0.25, 0.50, 1.00, 1.50, 2.00, 3.00]
W = [0.0, 0.25, 1.0]


def all_states(n, k=4):
    """Every configuration in {0..k-1}^n, in the same order the enum npz uses."""
    idx = np.arange(k ** n)
    S = np.empty((k ** n, n), dtype=np.int8)
    for r in range(n):                      # residue 0 is the most significant digit
        S[:, r] = (idx // (k ** (n - 1 - r))) % k
    return S


def hamil_over_enumeration(pdb, chunk=8192):
    """Native-free structural objective evaluated on the COMPLETE enumeration."""
    t = {x["pdb"]: x for x in I.targets()}[pdb]
    T = H.Terms(pdb, t["seq"], t["n"], t["fold"])
    S = all_states(T.n, T.k)
    ep = np.empty(len(S)); ed = np.empty(len(S))
    for a in range(0, len(S), chunk):
        b = S[a:a + chunk].astype(int)
        ep[a:a + chunk] = T.e_prior(b)
        ed[a:a + chunk] = T.e_disto(b, chunk=chunk)
    return H._z(ep), H._z(ed), S


def _verify_order(pdb, S, z):
    """The enum npz and our state enumeration must index the SAME configurations.

    Checked by rebuilding the snap configuration's index rather than trusting the ordering
    convention -- an off-by-one in digit order would silently pair every objective value
    with the wrong RMSD and produce a beautiful, meaningless result.
    """
    snap = np.asarray(z.snap_states, int)
    k, n = int(z.k), int(z.n)
    idx = 0
    for r in range(n):
        idx = idx * k + int(snap[r])
    ok = bool(np.array_equal(S[idx], snap))
    return ok, idx, int(z.snap_index)


def run(targets=ENUM, gaps=GAPS, ws=W):
    rows = {}
    for pdb in targets:
        z = E.enum(pdb)
        R = np.asarray(z.rmsd, float)
        ep, ed, S = hamil_over_enumeration(pdb)

        ok, idx, recorded = _verify_order(pdb, S, z)
        assert ok and idx == recorded, (
            f"{pdb}: state ordering does not match the enumeration cache "
            f"(rebuilt {idx}, recorded {recorded}) -- objective values would be paired "
            f"with the wrong RMSDs")

        objs = {"legacy": np.asarray(z.legacy, float),
                "prior_1local": np.asarray(z.prior, float)}
        for w in ws:
            objs[f"hamil_w{w:g}"] = (1.0 - w) * ep + w * ed

        rows[pdb] = {}
        for name, Ev in objs.items():
            rng = np.random.default_rng(0)
            rows[pdb][name] = {
                "pair_acc": {f"{g:g}": E.pair_accuracy(Ev, R, n_pairs=300000,
                                                       rng=np.random.default_rng(7),
                                                       min_gap=g) for g in gaps},
                "argmin_rmsd": E.argmin_rmsd(Ev, R),
                "decile_rho": E.decile_rho(Ev, R),
            }
        print(f"  {pdb} done")

    agg = {}
    for name in rows[targets[0]]:
        agg[name] = {
            "pair_acc": {f"{g:g}": float(np.mean([rows[p][name]["pair_acc"][f"{g:g}"]
                                                  for p in targets])) for g in gaps},
            "argmin_rmsd": float(np.mean([rows[p][name]["argmin_rmsd"] for p in targets])),
            "decile_rho": float(np.mean([rows[p][name]["decile_rho"] for p in targets])),
        }

    # the number that decides it: the gap at which each objective first reaches 55%
    thr = {}
    for name, a in agg.items():
        hit = None
        for g in gaps:
            if a["pair_acc"][f"{g:g}"] >= 0.55:
                hit = g
                break
        thr[name] = hit

    out = {"targets": targets, "gaps": gaps, "per_target": rows, "aggregate": agg,
           "gap_for_55pct": thr,
           "ener_reference": {"leg_steric": 0.58, "legacy": 1.31, "amber": 1.58,
                              "floor": "no physical objective exceeds 0.511 below 0.25 A"}}
    with open(os.path.join(RESULTS, "discrim.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_discrim", out, n_expected=len(targets))

    print(f"\nPairwise ranking accuracy vs quality gap, full enumeration, {len(targets)} "
          f"targets\n")
    hdr = "".join(f"{g:>9.2f}" for g in gaps)
    print(f"{'objective':<18}{hdr}{'55% at':>10}{'argmin':>9}")
    for name, a in agg.items():
        line = "".join(f"{a['pair_acc'][f'{g:g}']:>9.3f}" for g in gaps)
        t = thr[name]
        ts = f"{t:.2f} A" if t is not None else "never"
        print(f"{name:<18}{line}{ts:>10}{a['argmin_rmsd']:>9.3f}")
    print("\nENER reference on the physical objectives: leg_steric 0.58 A, Legacy 1.31 A, "
          "AMBER 1.58 A; no physical objective exceeds 0.511 below a 0.25 A gap.")
    return out


if __name__ == "__main__":
    run()
