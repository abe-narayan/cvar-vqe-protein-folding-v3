"""s16/verify_cosmag.py -- VERIFY workstream.  THE STATISTIC csteer.py QUOTES IS NOT THE ONE
ITS OWN LAW CONSUMES.

`s16/CSTEER_FINDINGS.md` section 2 reports the SIGNED mean cosine per direction and concludes
"from the best available starting structure, every native-free direction has a cosine with the
true residual that is statistically indistinguishable from zero".

But the law the same module derives, in its own docstring, is

    ||r - t v||^2 minimised at  t* = (v . r) = c ||r||     giving   ||r_new|| = ||r|| sqrt(1 - c^2)

`t*` carries the sign of `c`, and the payoff `sqrt(1 - c^2)` depends on `c^2`.  **The law is
two-sided.**  A direction with c = -0.4 is worth exactly as much as one with c = +0.4: you step
backwards along it.  So the quantity that has to be non-zero for the mechanism to pay is |c|,
and averaging the SIGNED c over targets destroys precisely the information the law uses.

Two consequences are tested here:

  C1  |c| against a matched random control, paired per target, fold-clustered.  |c| is
      positively biased for every direction, so the random control -- not zero -- is the
      reference.  E|cos| for an isotropic direction in 3n dimensions is sqrt(2/(pi*3n)) ~ 0.13
      at n = 13, which is what the control returns.

  C2  ZERO-INFORMATION REFERENCES on the same statistic (AUDIT's standing requirement): a
      constant ideal alpha-helix, beta-strand and polyproline-II, and an arbitrary retrieval
      window.  If one of those matches or beats the sprint's channels on |c|, section 2's
      framing dies the same way Sprint 15's 0.390 did.

  C3  THE LADDER IS ONE-SIDED.  `csteer.STEPS` runs 0.0 .. 1.0 with no negative rung, so no arm
      could ever exploit a negative cosine even where one exists.  With signed c ~ 0 and
      |c| ~ 0.3 this means roughly half the targets need a step the ladder cannot take.  The
      fraction is counted here.

Everything is reported as an ORACLE CEILING: the sign of c and the magnitude ||r|| are both
native-derived, so `base * (1 - sqrt(1 - c^2))` is what a perfect step selector would buy, never
a deployable effect.

Writes `s16/results/verify_cosmag.json`.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

HELIX = (np.deg2rad(-57.0), np.deg2rad(-47.0))
STRAND = (np.deg2rad(-139.0), np.deg2rad(135.0))
PPII = (np.deg2rad(-75.0), np.deg2rad(145.0))


def _sup(X, T):
    return I.superpose_batch(np.asarray(X, float)[None], np.asarray(T, float))[0]


def _cos(d, r):
    nd = float(np.linalg.norm(d)); nr = float(np.linalg.norm(r))
    if nd < 1e-12 or nr < 1e-12:
        return 0.0
    return float((d.ravel() @ r.ravel()) / (nd * nr))


def c1(out):
    d = json.load(open(os.path.join(RESULTS, "csteer.json")))
    rows = d["rows"]
    folds = np.asarray([r["fold"] for r in rows], int)
    pdbs = [r["pdb"] for r in rows]
    print("=" * 100)
    print("C1  THE ARTEFACT'S OWN COSINES, READ AS THE LAW READS THEM (|c|, not signed c)")
    print("=" * 100)
    res = {}
    for start in ("avg", "proj", "fit"):
        base = np.asarray([r["cells"][start]["base"] for r in rows], float)
        R = np.abs(np.asarray([[r["cells"][start]["arms"][f"rand{k}"]["cos"] for k in range(3)]
                               for r in rows], float)).mean(1)
        print(f"\n  START = {start}   mean {base.mean():.3f} A     "
              f"(random control |c| = {R.mean():.3f})")
        print(f"  {'direction':<14}{'signed c':>10}{'|c|':>8}"
              f"{'|c| - random, paired [95% CI]':>34}{'W/L':>9}"
              f"{'ORACLE two-sided ceiling':>26}")
        for dn in [k for k in rows[0]["cells"][start]["arms"] if not k.startswith("rand")]:
            c = np.asarray([r["cells"][start]["arms"][dn]["cos"] for r in rows], float)
            a = np.abs(c)
            ci = I.paired(a, R, folds=folds, names=pdbs)
            gain = base * (1.0 - np.sqrt(np.maximum(0.0, 1.0 - a ** 2)))
            w = int((a > R).sum()); l = int((a < R).sum())
            tag = "  ORACLE" if dn.startswith("ORACLE") else ""
            print(f"  {dn:<14}{c.mean():>+10.3f}{a.mean():>8.3f}"
                  f"{ci['mean_diff']:>+16.3f} [{ci['ci95'][0]:+.3f},{ci['ci95'][1]:+.3f}]"
                  f"{f'{w}/{l}':>9}{base.mean() - gain.mean():>16.3f} A{tag}")
            res.setdefault(start, {})[dn] = {
                "signed_mean": float(c.mean()), "abs_mean": float(a.mean()),
                "abs_minus_random": ci["mean_diff"], "ci95": ci["ci95"], "wl": [w, l],
                "oracle_two_sided_rmsd": float((base - gain).mean()),
                "frac_negative_c": float((c < 0).mean())}
        print(f"  requirement from the incumbent 3.204 A: |c| = 0.615 for 2.5 A, "
              f"0.781 for 2.0 A")
    out["c1"] = res

    print("\n" + "=" * 100)
    print("C3  THE LADDER IS ONE-SIDED: csteer.STEPS = 0.0 .. 1.0, no negative rung.")
    print("=" * 100)
    for start in ("avg", "proj", "fit"):
        line = f"  start {start:<6}"
        for dn in [k for k in rows[0]["cells"][start]["arms"]
                   if not k.startswith("rand") and not k.startswith("ORACLE")]:
            c = np.asarray([r["cells"][start]["arms"][dn]["cos"] for r in rows], float)
            line += f"  {dn} {float((c < 0).mean()):.2f}"
        print(line + "   <- fraction of targets whose optimal step is NEGATIVE")
    print("  On those targets no rung of the ladder can help, whatever the direction carries.")


def c2(out):
    from s14 import avgspace as AV
    from s15 import distcal as C
    print("\n" + "=" * 100)
    print("C2  ZERO-INFORMATION REFERENCES on |c|, from the best start (`avg`)")
    print("=" * 100)
    tg = I.targets()
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    pdbs = [t["pdb"] for t in tg]
    data = C.gather(tg)
    names = ["med", "helix", "strand", "ppii", "poolwin", "rand"]
    cos = {k: [] for k in names}
    base = []
    for t in tg:
        p = t["pdb"]; d = data[p]
        n = int(d["n"]); nat = np.asarray(d["nat"], float)
        W = np.asarray(AV.top75_windows(p)[0], float)
        P = I.pairwise_rmsd(W)
        avg_c, b_med = I.coordinate_average(W, P)
        S = np.asarray(avg_c, float)
        r = _sup(nat, S) - S
        base.append(float(I.ca_rmsd(S, nat)))
        rng = SD.stable_rng(p, "verify_csteer")
        dirs = {
            "med": np.asarray(W[b_med], float),
            "helix": I.build_ca(np.full(n, HELIX[0]), np.full(n, HELIX[1])),
            "strand": I.build_ca(np.full(n, STRAND[0]), np.full(n, STRAND[1])),
            "ppii": I.build_ca(np.full(n, PPII[0]), np.full(n, PPII[1])),
            "poolwin": np.asarray(W[int(rng.integers(0, len(W)))], float)}
        for nm, T in dirs.items():
            cos[nm].append(_cos(_sup(T, S) - S, r))
        z = rng.standard_normal(S.shape); z -= z.mean(0)
        cos["rand"].append(_cos(z, r))
    base = np.asarray(base, float)
    Rc = np.abs(np.asarray(cos["rand"], float))
    kinds = {"med": "native-free channel", "helix": "ZERO-INFORMATION",
             "strand": "ZERO-INFORMATION", "ppii": "ZERO-INFORMATION",
             "poolwin": "ZERO-INFO (arbitrary retrieval window)",
             "rand": "random control"}
    print(f"  n = {len(tg)}, start = `avg`, mean {base.mean():.3f} A")
    print(f"  {'direction':<10}{'kind':<40}{'signed c':>10}{'|c|':>8}"
          f"{'|c| - random [95% CI]':>30}{'ORACLE ceiling':>16}")
    res = {}
    for nm in names:
        c = np.asarray(cos[nm], float); a = np.abs(c)
        ci = I.paired(a, Rc, folds=folds, names=pdbs)
        gain = base * (1.0 - np.sqrt(np.maximum(0.0, 1.0 - a ** 2)))
        res[nm] = {"signed_mean": float(c.mean()), "abs_mean": float(a.mean()),
                   "abs_minus_random": ci["mean_diff"], "ci95": ci["ci95"],
                   "oracle_two_sided_rmsd": float((base - gain).mean()),
                   "per_target_abs": a.tolist()}
        print(f"  {nm:<10}{kinds[nm]:<40}{c.mean():>+10.3f}{a.mean():>8.3f}"
              f"{ci['mean_diff']:>+14.3f} [{ci['ci95'][0]:+.3f},{ci['ci95'][1]:+.3f}]"
              f"{(base - gain).mean():>14.3f} A")
    out["c2"] = res
    print("\n  A CONSTANT IDEAL BETA-STRAND -- which knows nothing about the target beyond its")
    print("  length -- carries MORE |c| from `avg` than any channel the sprint built.")


def main():
    out = {}
    c1(out)
    c2(out)
    with open(os.path.join(RESULTS, "verify_cosmag.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote s16/results/verify_cosmag.json")


if __name__ == "__main__":
    main()
