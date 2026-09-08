"""s16/diversity.py -- WHY the VQE ensembles read out worse, tested against a named identity.

`s16/integrate.py` produced a result its own diversity block does not explain by the
programme's standing law.  At m = 75 the six generators have essentially the SAME mean member
RMSD -- 3.886, 3.943, 3.903, 3.917, 3.897, 3.895, a spread of 0.057 A -- and yet their
coordinate-average readouts span 0.30 A (uniform 3.214 to cvar0.25 3.515).  The recorded law
"the terminal operator consumes the set MEAN, not the set BEST" therefore does not account
for it: the set means are equal.

THE HYPOTHESIS, and it is not ours.  The LIT workstream has just identified the programme's
"fusion law" as the **Krogh-Vedelsby (1995) ambiguity decomposition**, which for an ensemble
average says

    ||avg - native||^2  =  mean_i ||x_i - native||^2  -  mean_i ||x_i - avg||^2
    ensemble error       =  average error             -  DIVERSITY

If that is what is operating, then generators with equal average error must differ in
readout exactly by their diversity, and the VQE's readout penalty IS its concentration --
measured, not inferred.  Note the identity needs the QUADRATIC mean of the member errors;
LIT found the Sprint 15 record using the arithmetic mean, which is a live correction being
handled by the RETRACT workstream.  Both are printed here so the substitution's cost is
visible on this instrument too.

WHAT IS ACTUALLY IN QUESTION.  The identity is exact in any inner-product space.  It is NOT
obviously exact here, because every term is measured after Kabsch superposition, which is a
nonlinear operation, and because the coordinate average superposes members onto the ensemble
MEDOID rather than onto the native.  So the test has real content: if the residual is large,
the decomposition does not apply to this pipeline's operator and may not be used to explain
it.  That residual is the module's first output.

No AMBER and no Legacy are computed here -- this module asks a geometric question only, and
skipping them makes it cheap enough to run beside four agents.

ORACLE.  The native is read only to score.  Diversity itself is NATIVE-FREE, which is the
point: if diversity predicts the readout, a generator can be chosen without any native.
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
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I              # noqa: E402
from s14 import vqe_lib as V                 # noqa: E402
from s15 import seed as SD                   # noqa: E402
from s16.integrate import (ALPHAS, BUDGET, MS, SEEDS, SHOTS,   # noqa: E402
                           _ca_of, _ensembles, _pick)


def run(targets=None, seeds=SEEDS):
    from s13 import qarch_lib as QL
    tg = list(targets or V.ENUM_TARGETS)
    rows = []
    for ti, pdb in enumerate(tg):
        en = V.Enum(pdb)
        sp = QL.Space(pdb, en.k, seq=en.seq, n=en.n, fold=en.fold)
        nat = np.asarray(sp.nat_ca, float)
        for seed in seeds:
            for gname, seen in _ensembles(en, seed).items():
                for m in MS:
                    idx = _pick(en, seen, m)
                    W, _phi, _psi = _ca_of(en, idx)
                    P = I.pairwise_rmsd(W)
                    avg, b = I.coordinate_average(W, P)
                    #: every member is put in the SAME frame the average was built in, so the
                    #: decomposition's three terms are measured in one frame rather than three.
                    Wm = I.superpose_batch(W, W[b])
                    nn = float(len(nat))
                    per = np.array([I.ca_rmsd(w, nat) for w in W])
                    div2 = float(np.mean(((Wm - avg) ** 2).sum(axis=(1, 2)) / nn))
                    #: the member errors must be measured in the SAME frame as the diversity
                    #: for the identity to have a chance; also keep the free-superposition
                    #: version, which is what every other table in the sprint reports.
                    natm = I.superpose_batch(nat[None], avg)[0]
                    err2 = float(np.mean(((Wm - natm) ** 2).sum(axis=(1, 2)) / nn))
                    rows.append({
                        "pdb": pdb, "seed": int(seed), "gen": gname, "m": int(m),
                        "readout": float(I.ca_rmsd(avg, nat)),
                        "mean_member_arith": float(per.mean()),
                        "mean_member_quad": float(np.sqrt((per ** 2).mean())),
                        "err2_common_frame": err2, "div2": div2,
                        "kv_pred": float(np.sqrt(max(err2 - div2, 0.0))),
                        "diversity": float(np.sqrt(div2)),
                        "mean_pair_rmsd": float(P[np.triu_indices(len(W), 1)].mean()),
                    })
        print(f"  {ti+1}/{len(tg)} {pdb}", flush=True)
        json.dump({"rows": rows}, open(os.path.join(RESULTS, "diversity.json"), "w"))
    report(rows)
    return rows


def _boot(d, rng, B=4000):
    d = np.asarray(d, float); k = len(d)
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows):
    rng = SD.stable_rng("diversity", "report")
    gens = ["cvar0.1", "cvar0.25", "cvar0.5", "cvar1.0", "bestofN", "uniform"]
    print("\nA. DOES THE KROGH-VEDELSBY IDENTITY HOLD THROUGH THIS OPERATOR?")
    print("   (readout vs sqrt(mean member error^2 - diversity^2), all in the medoid frame)")
    print(f"  {'m':>4}{'readout':>10}{'KV prediction':>15}{'residual [95% CI]':>28}{'max |res|':>11}")
    for m in MS:
        sel = [r for r in rows if r["m"] == m]
        d = [r["readout"] - r["kv_pred"] for r in sel]
        a, lo, hi = _boot(d, rng)
        print(f"  {m:>4}{np.mean([r['readout'] for r in sel]):>10.3f}"
              f"{np.mean([r['kv_pred'] for r in sel]):>15.3f}"
              f"{a:>+16.4f} [{lo:+.4f},{hi:+.4f}]{max(abs(x) for x in d):>11.4f}")

    print("\n   the arithmetic-vs-quadratic mean substitution, on this instrument")
    for m in MS:
        sel = [r for r in rows if r["m"] == m]
        d = [r["mean_member_quad"] - r["mean_member_arith"] for r in sel]
        a, lo, hi = _boot(d, rng)
        print(f"  m = {m:<4} quadratic minus arithmetic mean member RMSD "
              f"{a:+.3f} [{lo:+.3f},{hi:+.3f}]  max {max(d):+.3f}")

    print("\nB. EQUAL AVERAGE ERROR, UNEQUAL DIVERSITY -- is diversity the whole story?")
    for m in MS:
        print(f"  --- m = {m}")
        print(f"  {'generator':<12}{'readout':>9}{'mean member':>13}{'diversity':>11}"
              f"{'mean pair RMSD':>16}")
        for g in gens:
            sel = [r for r in rows if r["gen"] == g and r["m"] == m]
            if not sel:
                continue
            print(f"  {g:<12}{np.mean([r['readout'] for r in sel]):>9.3f}"
                  f"{np.mean([r['mean_member_quad'] for r in sel]):>13.3f}"
                  f"{np.mean([r['diversity'] for r in sel]):>11.3f}"
                  f"{np.mean([r['mean_pair_rmsd'] for r in sel]):>16.3f}")

    print("\nC. IS DIVERSITY A NATIVE-FREE PREDICTOR OF THE READOUT?  (rank corr within m)")
    for m in MS:
        sel = [r for r in rows if r["m"] == m]
        x = np.array([r["mean_pair_rmsd"] for r in sel])
        y = np.array([r["readout"] for r in sel])
        rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
        rho = float(np.corrcoef(rx, ry)[0, 1])
        #: WITHIN TARGET, which is the only comparison a generator choice ever faces --
        #: across targets this correlation is dominated by target difficulty.
        wt = []
        for p in sorted({r["pdb"] for r in sel}):
            s2 = [r for r in sel if r["pdb"] == p]
            if len(s2) < 4:
                continue
            a = np.array([r["mean_pair_rmsd"] for r in s2])
            b = np.array([r["readout"] for r in s2])
            wt.append(float(np.corrcoef(np.argsort(np.argsort(a)),
                                        np.argsort(np.argsort(b)))[0, 1]))
        print(f"  m = {m:<4} pooled rho {rho:+.3f}   within-target rho "
              f"{np.mean(wt):+.3f} (median {np.median(wt):+.3f}, {len(wt)} targets)")
    print("\n  A NEGATIVE rho means more diverse ensembles read out BETTER, which is what the"
          "\n  ambiguity decomposition predicts and what would make diversity a native-free"
          "\n  generator-selection criterion.")


if __name__ == "__main__":
    run()
