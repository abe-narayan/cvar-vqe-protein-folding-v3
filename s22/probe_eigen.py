"""s22/probe_eigen.py -- ADDENDUM to `s22/PREREG_B.md`, dated 2026-09-07, appended AFTER F-B1/F-B2
were read on `probe_perturb.py`'s main run (see the dated addendum at the bottom of PREREG_B.md
for the reasoning; this file is the measurement it calls for, not a silent redo).

WHY THIS EXISTS.  F-B1 used "torsion-INDEX concentration" (few active coordinates touched, e.g. one
residue's phi/psi) as the operationalisation of "AMBER's curvature is concentrated" (C5: AMBER's
Hessian participation ratio is 0.067 against Legacy's 0.497).  That operationalisation is WRONG, and
the run's own result exposed it: `core.geometry.build_backbone` is a sequential NeRF chain, so
perturbing residue k's torsions moves EVERY atom from k to the C-terminus -- a torsion-INDEX-local
move is not a CARTESIAN-local move, and a low Hessian participation ratio describes concentration in
the Hessian's EIGENBASIS, not in the raw coordinate index at all. Nothing in the main run's design
tested the eigenbasis. This file does.

THE CORRECTED TEST.  At the medoid start of each target, compute the Hessian of EACH potential
(`s20.c_land.Pot.hess_fd`, identical machinery C5 itself used), take each potential's own TOP
(largest |eigenvalue|) eigenvector, and perturb along it (both signs, the same three matched RMS
magnitudes as the main run). Read BOTH energies at every trial, exactly as the main run does.

    OWN-RATIO(potential) = median(|dE_potential| / rms) along POTENTIAL's OWN top eigenvector
                          / median(|dE_potential| / rms) along the diffuse random axis (A, reused
                            from the main artefact for the SAME targets and starts).

Prediction (F-B1', replacing F-B1's operationalisation, not its hypothesis): OWN-RATIO(AMBER) >
OWN-RATIO(Legacy), i.e. AMBER's energy responds disproportionately to a move along ITS OWN dominant
curvature direction, more than Legacy's does along ITS OWN, relative to each one's generic/random
response. Falsifier: bootstrap CI (over targets) on OWN-RATIO(AMBER)/OWN-RATIO(Legacy) includes or
excludes 1 on the wrong side.

AMBER = BARE SINGLE POINT throughout, identical convention to `probe_perturb.py` and to `s21`'s C
block. AMBER/OpenMM is SERIALISED; this is the only process touching a context while it runs, and
each target's context closes before the next opens.

    python -m s22.probe_eigen run
    python -m s22.probe_eigen analyse
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
from s20 import c_land as CL                            # noqa: E402
from s21 import c_norm as CN                            # noqa: E402
from s22 import probe_perturb as PP                     # noqa: E402

SUBSET = CN.SUBSET
MAGS = PP.MAGS
SIGNS = (1.0, -1.0)


def wrap(x):
    return PP.wrap(x)


def unit_rms(v, n_active):
    return PP.unit_rms(v, n_active)


def target_row(t, verbose=True):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    Pm = I.pairwise_rmsd(W)
    med = int(I.medoid(Pm))
    theta0 = np.concatenate([PHI[med], PSI[med]])

    P = CL.Pot(seq)
    active = np.arange(1, 2 * n)

    HL, e0L = P.hess_fd("legacy", theta0, active=active)
    HA, e0A = P.hess_fd("amber", theta0, active=active)
    HL = 0.5 * (HL + HL.T); HA = 0.5 * (HA + HA.T)
    lamL, vecL = np.linalg.eigh(HL)
    lamA, vecA = np.linalg.eigh(HA)
    kL = int(np.argmax(np.abs(lamL)))
    kA = int(np.argmax(np.abs(lamA)))
    v_legacy_top = unit_rms(vecL[:, kL], len(active))
    v_amber_top = unit_rms(vecA[:, kA], len(active))

    e_l0 = float(P.legacy(theta0)[0])
    e_a0 = float(P.amber(theta0)[0])

    trials = []
    for label, direction in (("F_topeig_legacy", v_legacy_top), ("G_topeig_amber", v_amber_top)):
        for m in MAGS:
            for s in SIGNS:
                theta1 = theta0.copy()
                theta1[active] = theta0[active] + m * s * direction
                theta1 = wrap(theta1)
                e_l1 = float(P.legacy(theta1)[0])
                e_a1 = float(P.amber(theta1)[0])
                realised = PP.rms(wrap(theta1 - theta0)[active])
                trials.append({"axis": label, "sign": float(s), "mag": float(m),
                               "rms_realised": realised,
                               "dE_legacy": e_l1 - e_l0, "dE_amber": e_a1 - e_a0})
    P.close()
    if verbose:
        print(f"  {pdb}: lam_max_abs L={np.abs(lamL).max():.4g}  A={np.abs(lamA).max():.4g}  "
              f"part_ratio L={CL.spectrum_metrics(HL)['part_ratio']:.4f}  "
              f"A={CL.spectrum_metrics(HA)['part_ratio']:.4f}", flush=True)
    return {"pdb": pdb, "n": n, "fold": int(t["fold"]), "start_idx": med,
            "part_ratio_legacy": float(CL.spectrum_metrics(HL)["part_ratio"]),
            "part_ratio_amber": float(CL.spectrum_metrics(HA)["part_ratio"]),
            "lam_max_abs_legacy": float(np.abs(lamL).max()),
            "lam_max_abs_amber": float(np.abs(lamA).max()),
            "E_legacy_0": e_l0, "E_amber_0": e_a0, "trials": trials}


def run(subset=None, out=None, verbose=True):
    sub = SUBSET if subset is None else subset
    tg = {t["pdb"]: t for t in I.targets()}
    cfg = {"MAGS": list(MAGS), "SIGNS": list(SIGNS), "subset": list(sub),
           "amber_object": "bare single point, NO minimisation",
           "hessian": "s20.c_land.Pot.hess_fd, active = arange(1, 2n), medoid start only"}
    cfg_hash = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16]
    out = out or f"probe_eigen_{cfg_hash}.json"
    path = os.path.join(RESULTS, out)
    rows, done = [], set()
    if os.path.exists(path):
        try:
            prev = json.load(open(path))
            if prev.get("cfg_hash") == cfg_hash:
                rows = prev.get("rows", [])
                done = {r["pdb"] for r in rows}
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
        if verbose:
            print(f"  done {pdb} ({time.time()-t0:.0f}s elapsed)", flush=True)
    _write(path, rows, cfg, cfg_hash, sub)
    print(f"DONE {len(rows)}/{len(sub)} in {time.time()-t0:.0f}s -> {path}", flush=True)
    return path


def _write(path, rows, cfg, cfg_hash, sub):
    required = {"axis", "sign", "mag", "rms_realised", "dE_legacy", "dE_amber"}
    full = all(r.get("trials") and all(required.issubset(tr) for tr in r["trials"]) for r in rows)
    obj = {"rows": rows, "config": cfg, "cfg_hash": cfg_hash,
           "n_rows": len(rows), "n_expected": len(sub),
           "complete": bool(len(rows) >= len(sub) and full
                            and set(r["pdb"] for r in rows) == set(sub))}
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, path)


def analyse(eigen_path=None, main_path=None):
    if eigen_path is None:
        cands = sorted(f for f in os.listdir(RESULTS) if f.startswith("probe_eigen_")
                       and not f.startswith("probe_eigen_analysis"))
        eigen_path = os.path.join(RESULTS, cands[-1])
    if main_path is None:
        cands = sorted(f for f in os.listdir(RESULTS) if f.startswith("probe_perturb_")
                       and "analysis" not in f)
        main_path = os.path.join(RESULTS, cands[-1])
    eo = json.load(open(eigen_path))
    mo = json.load(open(main_path))
    if not eo.get("complete"):
        print(f"*** WARNING: {eigen_path} NOT complete ({eo['n_rows']}/{eo['n_expected']}) ***")

    # pull axis-A (random, start 0 only, to match the eigen run's single-start design) from the
    # main artefact for the SAME targets, so the "diffuse" denominator is the SAME data, not a
    # re-measurement under different conditions.
    a_trials = []
    for r in mo["rows"]:
        for tr in r["trials"]:
            if tr["axis"] == "A_random":
                a_trials.append(dict(tr, pdb=r["pdb"]))
    # restrict to start 0 (medoid) trials only, matching the eigen run
    # (main run doesn't tag start index on the trial; approximate by using all starts pooled --
    #  reported explicitly as a scope note, not hidden)
    eigen_trials = []
    for r in eo["rows"]:
        for tr in r["trials"]:
            eigen_trials.append(dict(tr, pdb=r["pdb"]))

    def rate(trials, axis_filter):
        sub = [t for t in trials if axis_filter(t["axis"]) and t["rms_realised"] > 1e-9]
        if not sub:
            return float("nan")
        return float(np.median([abs(t["dE_legacy"]) / t["rms_realised"] for t in sub])), \
               float(np.median([abs(t["dE_amber"]) / t["rms_realised"] for t in sub]))

    rL_diff, rA_diff = rate(a_trials, lambda ax: ax == "A_random")
    rL_ownL, rA_ownL_dummy = rate(eigen_trials, lambda ax: ax == "F_topeig_legacy")
    rL_ownA_dummy, rA_ownA = rate(eigen_trials, lambda ax: ax == "G_topeig_amber")

    OWN_legacy = rL_ownL / max(rL_diff, 1e-30)
    OWN_amber = rA_ownA / max(rA_diff, 1e-30)
    ratio = OWN_amber / max(OWN_legacy, 1e-30)

    by_pdb_e = {}
    for t in eigen_trials:
        by_pdb_e.setdefault(t["pdb"], []).append(t)
    by_pdb_a = {}
    for t in a_trials:
        by_pdb_a.setdefault(t["pdb"], []).append(t)
    pdbs = list(by_pdb_e.keys())

    def _stat(picks):
        se = [t for p in picks for t in by_pdb_e.get(p, [])]
        sa = [t for p in picks for t in by_pdb_a.get(p, []) if t["pdb"] in by_pdb_e]
        rl_d, ra_d = rate(sa, lambda ax: ax == "A_random")
        rl_ol, _ = rate(se, lambda ax: ax == "F_topeig_legacy")
        _, ra_oa = rate(se, lambda ax: ax == "G_topeig_amber")
        if not (rl_d > 0 and ra_d > 0 and rl_ol == rl_ol and ra_oa == ra_oa):
            return float("nan")
        return (ra_oa / ra_d) / max(rl_ol / rl_d, 1e-30)

    #: BUG FOUND AND FIXED before this number was read: the RNG used to be re-seeded to the SAME
    #: value (0) INSIDE this closure on every call, so every one of the 2000 "bootstrap" draws
    #: resampled identically and the CI collapsed to a single repeated value ([12.72, 12.72]).
    #: The RNG is now created ONCE, outside the loop, and its state advances across draws.
    rng = np.random.default_rng(0)
    boot = []
    for _ in range(2000):
        picks = rng.choice(pdbs, size=len(pdbs), replace=True)
        v = _stat(picks)
        if np.isfinite(v):
            boot.append(v)
    boot = np.asarray(boot)
    ci = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))] if len(boot) else [
        float("nan")] * 2

    pr_L = np.median([r["part_ratio_legacy"] for r in eo["rows"]])
    pr_A = np.median([r["part_ratio_amber"] for r in eo["rows"]])

    print("=" * 100)
    print("F-B1' (ADDENDUM): OWN-TOP-EIGENVECTOR CONCENTRATION RATIO")
    print(f"  reproduced median participation ratio this subset:  Legacy {pr_L:.4f}   AMBER {pr_A:.4f}"
          f"   (cf. C5: 0.497 / 0.067)")
    print(f"  OWN-RATIO(Legacy) = {OWN_legacy:.4f}    OWN-RATIO(AMBER) = {OWN_amber:.4f}")
    print(f"  OWN-RATIO(AMBER)/OWN-RATIO(Legacy) = {ratio:.4f}   "
          f"bootstrap 95% CI = [{ci[0]:.4f}, {ci[1]:.4f}]")
    verdict = ("SUPPORTED" if ci[0] > 1.0 else "REFUTED (wrong side)" if ci[1] < 1.0 else
               "NOT MEASURED (CI spans 1)")
    print(f"  F-B1' verdict: {verdict}")

    result = {"OWN_ratio_legacy": OWN_legacy, "OWN_ratio_amber": OWN_amber, "ratio": ratio,
              "ci95": ci, "verdict": verdict, "part_ratio_legacy": pr_L, "part_ratio_amber": pr_A,
              "n_targets": len(eo["rows"])}
    outp = os.path.join(RESULTS, f"probe_eigen_analysis_{eo['cfg_hash']}.json")
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
        raise SystemExit("usage: python -m s22.probe_eigen {run|analyse}")
