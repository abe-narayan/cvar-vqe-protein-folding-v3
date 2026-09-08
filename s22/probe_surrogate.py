"""s22/probe_surrogate.py -- PRIORITY 2: IS THERE A DEFENSIBLE AMBER SURROGATE?

Candidate: `s21.c_cont.MixPot.amber_bonded` -- the BONDED force-group subset (bond + angle +
torsion) of the SAME genuine ff14SB System, nonbonded and solvation switched off. It is a subset of
the real Hamiltonian, not a learned or fitted model, and `s21/c_cont.py`'s own docstring is explicit
that it is "NEVER called amber anywhere" (PREREG section 6 there). This file inherits that
discipline: the candidate is called `H_surrogate` or `amber_bonded` throughout, never `amber`, and
the genuine bare-single-point `H_AMBER` is computed on every trial and reported alongside, never
displaced (Pillar 3).

QUESTION.  `core/amber.py`'s own docstring measures `CustomGBForce` (nonbonded solvation) at ~98% of
one AMBER evaluation's cost. If the bonded subset alone tracks the FULL single point's RESPONSE to a
torsion perturbation well enough to preserve RANKING (not magnitude -- it structurally cannot: it is
missing 98% of the energy), it is a defensible cheap proxy for screening; if not, no surrogate
survives and that is reported as plainly as a positive would be.

Reuses the SAME 30-target subset and the SAME axis-A (diffuse) and axis-B (compactness) trial
constructions as `probe_perturb.py`, at the medoid start only, to keep this a bounded add-on rather
than a second full run. AMBER/OpenMM SERIALISED -- sole context holder while this runs.

    python -m s22.probe_surrogate run
    python -m s22.probe_surrogate analyse
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                       # noqa: E402
from s14.avgspace import top75_windows                # noqa: E402
from s15 import seed as SD                             # noqa: E402
from s21 import c_norm as CN                            # noqa: E402
from s21 import c_cont as CC                            # noqa: E402
from s22 import probe_perturb as PP                     # noqa: E402

SUBSET = CN.SUBSET
MAGS = PP.MAGS
N_REPEAT_RANDOM = PP.N_REPEAT_RANDOM
SIGNS = PP.SIGNS


def target_row(t, verbose=True):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    med = int(I.medoid(I.pairwise_rmsd(W)))
    theta0 = np.concatenate([PHI[med], PSI[med]])

    P = CC.MixPot(seq)
    from s22.probe_perturb import TargetGeom, build_directions, random_dir, wrap
    TG = TargetGeom(seq)
    active = np.arange(1, 2 * n)
    dirs = build_directions(theta0, n, TG, active)

    e_a0 = float(P.amber(theta0)[0])
    e_b0 = float(P.amber_bonded(theta0)[0])

    trials = []
    rngA = SD.stable_rng(pdb, "s22surrogateA")
    for m in MAGS:
        for r in range(N_REPEAT_RANDOM):
            d = random_dir(len(active), rngA)
            theta1 = wrap(theta0.copy())
            theta1[active] = theta0[active] + m * d
            theta1 = wrap(theta1)
            e_a1 = float(P.amber(theta1)[0])
            e_b1 = float(P.amber_bonded(theta1)[0])
            trials.append({"axis": "A_random", "mag": float(m), "repeat": int(r),
                           "dE_amber": e_a1 - e_a0, "dE_bonded": e_b1 - e_b0})
    for m in MAGS:
        for s in SIGNS:
            d = dirs["ext_u"] if s > 0 else dirs["hel_u"]
            theta1 = theta0.copy()
            theta1[active] = theta0[active] + m * s * d
            theta1 = wrap(theta1)
            e_a1 = float(P.amber(theta1)[0])
            e_b1 = float(P.amber_bonded(theta1)[0])
            trials.append({"axis": "B_compactness", "mag": float(m), "sign": float(s),
                           "dE_amber": e_a1 - e_a0, "dE_bonded": e_b1 - e_b0})
    P.close()
    if verbose:
        print(f"  {pdb}: E_amber_0={e_a0:.4g}  E_bonded_0={e_b0:.4g}  "
              f"bonded_frac_of_magnitude={abs(e_b0)/max(abs(e_a0),1e-9):.4f}", flush=True)
    return {"pdb": pdb, "n": n, "E_amber_0": e_a0, "E_bonded_0": e_b0, "trials": trials}


def run(subset=None, out=None, verbose=True):
    sub = SUBSET if subset is None else subset
    tg = {t["pdb"]: t for t in I.targets()}
    cfg = {"MAGS": list(MAGS), "N_REPEAT_RANDOM": N_REPEAT_RANDOM, "SIGNS": list(SIGNS),
           "subset": list(sub), "start": "medoid only",
           "surrogate": "s21.c_cont.MixPot.amber_bonded -- bond+angle+torsion force groups of the "
                        "SAME genuine ff14SB System, nonbonded/solvation OFF. NEVER called amber."}
    cfg_hash = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16]
    out = out or f"probe_surrogate_{cfg_hash}.json"
    path = os.path.join(RESULTS, out)
    rows, done = [], set()
    if os.path.exists(path):
        try:
            prev = json.load(open(path))
            if prev.get("cfg_hash") == cfg_hash:
                rows = prev.get("rows", []); done = {r["pdb"] for r in rows}
        except Exception:
            rows = []
    t0 = time.time()
    for pdb in sub:
        if pdb in done:
            continue
        if verbose:
            print(f"[{len(rows)+1}/{len(sub)}] {pdb} ...", flush=True)
        rows.append(CN.with_mem_retry(target_row, tg[pdb], verbose=verbose))
        _write(path, rows, cfg, cfg_hash, sub)
    _write(path, rows, cfg, cfg_hash, sub)
    print(f"DONE {len(rows)}/{len(sub)} in {time.time()-t0:.0f}s -> {path}", flush=True)
    return path


def _write(path, rows, cfg, cfg_hash, sub):
    required = {"axis", "mag", "dE_amber", "dE_bonded"}
    full = all(r.get("trials") and all(required.issubset(tr) for tr in r["trials"]) for r in rows)
    obj = {"rows": rows, "config": cfg, "cfg_hash": cfg_hash,
           "n_rows": len(rows), "n_expected": len(sub),
           "complete": bool(len(rows) >= len(sub) and full
                            and set(r["pdb"] for r in rows) == set(sub))}
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, path)


def analyse(path=None):
    if path is None:
        cands = sorted(f for f in os.listdir(RESULTS) if f.startswith("probe_surrogate_")
                       and "analysis" not in f)
        path = os.path.join(RESULTS, cands[-1])
    o = json.load(open(path))
    if not o.get("complete"):
        print(f"*** WARNING: {path} NOT complete ({o['n_rows']}/{o['n_expected']}) ***")
    from scipy.stats import spearmanr
    all_a, all_b = [], []
    frac = []
    per_axis = {}
    for r in o["rows"]:
        frac.append(abs(r["E_bonded_0"]) / max(abs(r["E_amber_0"]), 1e-9))
        for tr in r["trials"]:
            all_a.append(tr["dE_amber"]); all_b.append(tr["dE_bonded"])
            per_axis.setdefault(tr["axis"], {"a": [], "b": []})
            per_axis[tr["axis"]]["a"].append(tr["dE_amber"])
            per_axis[tr["axis"]]["b"].append(tr["dE_bonded"])
    rho_all = spearmanr(all_a, all_b).correlation
    print("=" * 100)
    print("PRIORITY 2 -- H_surrogate (amber_bonded) vs genuine bare-single-point H_AMBER")
    print(f"  median |E_bonded_0| / |E_amber_0| at baseline (magnitude the surrogate structurally "
          f"CANNOT capture): {np.median(frac):.4f}")
    print(f"  Spearman(dE_amber, dE_bonded) pooled over all trials n={len(all_a)}: {rho_all:.4f}")
    for ax, d in per_axis.items():
        rho = spearmanr(d["a"], d["b"]).correlation
        print(f"    axis {ax:<16} n={len(d['a']):>4}  Spearman = {rho:.4f}")
    verdict = ("a defensible RANKING surrogate for torsion-perturbation response (rho > 0.7)"
               if rho_all > 0.7 else
               "NOT a defensible surrogate even for ranking (rho <= 0.7)")
    print(f"\nVERDICT: {verdict}")
    result = {"median_bonded_fraction_of_magnitude": float(np.median(frac)),
              "spearman_pooled": float(rho_all),
              "spearman_per_axis": {ax: float(spearmanr(d["a"], d["b"]).correlation)
                                    for ax, d in per_axis.items()},
              "verdict": verdict, "n_targets": len(o["rows"])}
    outp = os.path.join(RESULTS, f"probe_surrogate_analysis_{o['cfg_hash']}.json")
    with open(outp, "w") as fh:
        json.dump(result, fh, indent=1)
    print(f"\nwritten {outp}")
    return result


if __name__ == "__main__":
    a = sys.argv[1:]
    cmd = a[0] if a else "run"
    if cmd == "run":
        run()
    elif cmd == "analyse":
        analyse()
    else:
        raise SystemExit("usage: python -m s22.probe_surrogate {run|analyse}")
