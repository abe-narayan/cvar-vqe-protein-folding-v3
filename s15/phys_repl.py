"""SPRINT 15, PHYS -- is the programme's ONE positive physics result above the noise floor?

THE RESULT UNDER TEST.  Sprint 14 (`s14/ener_avgrefine.py`, `s14/results/ener_avgrefine.json`)
measured that restrained AMBER ff14SB/GBn2 relaxation of the ALL-ATOM COORDINATE AVERAGE at
`k_restraint = 30` emits a structure **-0.022 A [-0.036, -0.009]** better than the incumbent
ideal-geometry projection of the same coordinate average, at valid geometry and with a higher
Ramachandran-favoured fraction.  It is the only positive physics result the programme has.

THE CHALLENGE.  The RESTRAINT workstream demonstrated an empirical FALSE-POSITIVE FLOOR of
~0.08 A for paired comparisons in this sprint's machinery, by running a comparison that is
provably zero (maximum likelihood against a constant-width Gaussian IS least squares) and
getting `+0.081 [+0.014, +0.169]` on one multi-start draw and `-0.003` on another.  -0.022 is a
quarter of that floor.

THE FIRST THING THIS MODULE ESTABLISHES, AND IT CHANGES THE QUESTION.
**The k = 30 pipeline contains no random number generator at all.**

    top75_windows      cache read, deterministic
    coordinate_average medoid + superposition, deterministic
    build_backbone_*   deterministic
    I.project          `core.project.STARTS` is a FIXED 4-tuple of canonical (phi, psi)
                       constants -- (-120,130), (-57,-47), (-139,135), (-75,145).
                       `multi=True` means "all four fixed starts", NOT "random starts".
    refine_coords      OpenMM LocalEnergyMinimizer, steps=0, tolerance 1.0, threads=1

So the mechanism that produced the 0.08 A floor -- a per-process-salted `hash()` seeding a
random multi-start draw -- **cannot act here**, and the 0.08 A figure may not simply be imported.
The floor has to be built for THIS pipeline.  `--mode determinism` verifies the claim empirically
by re-running targets in a fresh interpreter and comparing bit patterns.

THE THREE DRAW AXES, and why each is the analogue it claims to be.

  `canon`  draw 0.  Exactly the Sprint 14 construction.  Reproduces -0.022 or the module is
           measuring something else.

  `frame`  **THE EXACT NULL.**  Apply a random proper rotation and translation to the averaged
           backbone before handing it to BOTH arms.  Every quantity in the comparison is
           rigid-invariant BY CONSTRUCTION: ff14SB/GBn2 is translation- and rotation-invariant,
           the positional restraint is to the reference in the same frame, `fit_prior` minimises
           a Kabsch RMSD, and `ca_rmsd` is a Kabsch RMSD.  **The frame-to-frame difference is
           mathematically zero**, exactly as ML-against-a-constant-width-Gaussian is
           mathematically zero.  Anything it measures is the pipeline's own numerical floor.
           This is not hypothetical: `verify/amber_audit.json` records the MINIMISER moving
           0.148 kcal/mol under translation while the single point moves 1.3e-6.

  `startperm`
           **THE SECOND EXACT NULL, and the one that reproduces the 0.08 A MECHANISM.**  The
           `frame` null perturbs floating point only; both arms stay in the same basin, so it
           certifies the numerical floor and nothing more.  R1.5's floor came from something
           stronger -- a multi-start optimiser landing in DIFFERENT basins on different draws.
           The only multi-start in this pipeline is `fit_multi`, which loops over the four
           fixed `STARTS` and keeps the lowest objective with a STRICT `<`.  Permuting the
           order of those four starts is **zero by construction** (an argmin over an unordered
           set) -- unless two starts are within float error of each other, in which case the
           strict `<` breaks the tie BY ORDER and the emitted structure flips basins.  That is
           this project's own recorded trap ("an argmin over a tied signal read the ORACLE sort
           order"), turned into a null.  Arm B is untouched, so only arm A is recomputed.

  `boot`   **THE REPLICATION AXIS.**  Bootstrap the 75 windows that enter the average (75 draws
           with replacement, ~63% distinct).  Both arms consume the identical resampled set, so
           the comparison stays start-matched, but the object being projected and relaxed is a
           different, equally legitimate draw of the pipeline's input.  This is the closest
           available analogue of "an independent start draw": it perturbs what the two local
           optimisers are handed without privileging either arm.

A DEFECT IN THE SPRINT 14 REPORT THAT THIS MODULE FIXES.  `s14/ener_avgrefine.run_target` scores
arm A's RMSD on `pr["fit_ca"]` -- the **lam = 0** projection -- and arm A's Ramachandran fraction
on `pr["phi"], pr["psi"]` -- the **lam = 0.3 ramah-penalised** projection.  Those are two
different structures, and the penalised one is the arm that was explicitly optimised for
Ramachandran plausibility.  The validity half of the claim was therefore comparing AMBER against
a structure whose RMSD is not the one quoted.  `project_both` returns both arms' torsions so the
two statistics can be reported on the SAME structure.  Arm A's clash count was also reported as a
hard-coded 0.0 placeholder; it is measured here.

    python -m s15.phys_repl --mode determinism
    python -m s15.phys_repl --mode canon
    python -m s15.phys_repl --mode frame --draws 1,2,3
    python -m s15.phys_repl --mode boot  --draws 1,2,3,4,5
    python -m s15.phys_repl --report
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import json
import math
import sys
import time

import numpy as np

from s12 import instrument as I
from s14 import ener_lib as E
from s14.avgspace import top75_windows
from s14.ener_avgrefine import ATOMS, IDEAL, _kabsch, concentration_null, geom_of, rama_ok, torsions_of
from s15.seed import stable_rng

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(RESULTS, exist_ok=True)
K = 30.0

#: The pristine generic start set, captured once at import.  `startperm` permutes THIS,
#: never the live tuple, so every draw is a permutation of draw 0 rather than of its
#: predecessor.
from core import project as _pj  # noqa: E402
_STARTS0 = tuple(tuple(s) for s in _pj.STARTS)


# ------------------------------------------------------------------ construction
def averaged_backbone_from(W, PHI, PSI):
    """`s14.ener_avgrefine.averaged_backbone` refactored to take explicit window arrays.

    Identical arithmetic; the only change is that the window set is a parameter, so a
    bootstrap resample can be fed through the exact same construction.
    """
    from core import geometry as geo
    W = np.asarray(W, float)
    C_ca, b = I.coordinate_average(W)
    Wm = I.superpose_batch(W, W[b])
    full = geo.build_backbone_batch(PHI, PSI)
    bca = np.asarray(full["CA"], float)
    R, Am, Bm = _kabsch(bca, Wm)
    avg = {}
    for a in ATOMS:
        if a not in full:
            continue
        X = np.asarray(full[a], float)
        avg[a] = (np.einsum("bij,bnj->bni", R, X - Am) + Bm).mean(0)
    dev = float(np.abs(avg["CA"] - C_ca).max())
    return avg, C_ca, dev


def project_both(C, seq, fold, lam=0.3, maxiter=300):
    """`I.project`, but returning the lam = 0 arm's TORSIONS as well as its CA trace.

    `I.project` throws away `fit`'s torsions, which is why Sprint 14 had to score arm A's
    Ramachandran fraction on the penalised arm while quoting the unpenalised arm's RMSD.
    """
    from core import project as pj
    pen = pj.make_penalty("ramah", seq, int(fold))
    path = pj.lam_path(np.asarray(C, float), pen, (0.0, lam), maxiter=maxiter,
                       multi=True, grad="exact")
    fit, arm = path[0.0], path[lam]
    return dict(fit_ca=np.asarray(fit[0], float),
                fit_phi=np.asarray(fit[1], float), fit_psi=np.asarray(fit[2], float),
                ca=np.asarray(arm[0], float),
                phi=np.asarray(arm[1], float), psi=np.asarray(arm[2], float))


def random_rigid(rng):
    """A proper rotation (det = +1) and a translation, for the EXACT NULL."""
    Q, R = np.linalg.qr(rng.normal(size=(3, 3)))
    Q = Q * np.sign(np.diag(R))          # fix QR's sign ambiguity -> Haar on O(3)
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1.0                  # -> SO(3); a reflection is NOT a symmetry of ff14SB
    return Q, rng.normal(scale=10.0, size=3)


# ------------------------------------------------------------------ one target, one draw
def run_target(t, mode, draw, k=K, do_proj=True, do_amber=True):
    """One target, one draw.

    `do_proj=False` skips arm A.  Legitimate ONLY in `frame` mode, and only because arm A is
    frame-invariant BY CONSTRUCTION -- `fit_prior` minimises `kabsch_rmsd(build(phi,psi), C)`,
    which is a function of C's internal geometry alone, so a rigid motion of C cannot change
    the L-BFGS trajectory except through float error in the Kabsch step itself.  That claim is
    not assumed: `frame` draw 1 runs BOTH arms on all 126 targets and reports arm A's own
    frame-to-frame deviation, and only draws 2+ skip it.
    """
    from core import amber as am
    from core import geometry as geo
    import torsion_lib2 as tl2
    pdb, seq, fold = t["pdb"], t["seq"], int(t["fold"])
    W, PHI, PSI, u = top75_windows(pdb)
    m = len(W)
    rng = stable_rng("phys15", mode, int(draw), pdb)

    boot_idx = None
    perm = None
    if mode == "startperm":
        from core import project as pj
        perm = rng.permutation(len(_STARTS0))
        # permute the PRISTINE tuple, never the current one -- permuting in place would
        # compound across targets and the draw would stop being a permutation of draw 0.
        pj.STARTS = tuple(_STARTS0[i] for i in perm)
        rec_perm = [int(x) for x in perm]
    if mode == "boot":
        boot_idx = rng.integers(0, m, m)
        W, PHI, PSI = np.asarray(W)[boot_idx], np.asarray(PHI)[boot_idx], np.asarray(PSI)[boot_idx]

    avg, C_ca, dev = averaged_backbone_from(W, PHI, PSI)

    if mode == "frame":
        Q, tv = random_rigid(rng)
        avg = {a: v @ Q.T + tv for a, v in avg.items()}
        C_ca = C_ca @ Q.T + tv

    nat = u["nat_ca"]
    rec = dict(pdb=pdb, fold=fold, n=len(seq), mode=mode, draw=int(draw),
               n_windows=int(m), n_distinct=int(len(np.unique(boot_idx))) if boot_idx is not None else int(m),
               ca_identity_max_dev=dev,
               raw_avg_rmsd=float(I.ca_rmsd(avg["CA"], nat)),
               raw_avg_geom=geom_of(avg, seq))
    ph, ps = torsions_of(avg)
    rec["raw_avg_rama_ok"] = rama_ok(ph, ps)

    # ---- arm A, the incumbent projection.  BOTH lambda arms, scored consistently.
    if do_proj:
        t0p = time.time()
        pr = project_both(C_ca, seq, fold)
        rec["projected_rmsd"] = float(I.ca_rmsd(pr["fit_ca"], nat))       # the quoted arm (lam=0)
        rec["projected_lam_rmsd"] = float(I.ca_rmsd(pr["ca"], nat))       # the ramah-penalised arm
        rec["projected_rama_ok"] = rama_ok(pr["fit_phi"], pr["fit_psi"])  # MATCHED to the quoted arm
        rec["projected_lam_rama_ok"] = rama_ok(pr["phi"], pr["psi"])      # what S14 reported instead
        bbA = geo.build_backbone(pr["fit_phi"], pr["fit_psi"])
        rec["projected_geom"] = geom_of({a: np.asarray(bbA[a], float) for a in ATOMS if a in bbA}, seq)
        rec["projected_wall"] = round(time.time() - t0p, 2)

    if perm is not None:
        rec["start_perm"] = rec_perm
        _pj.STARTS = _STARTS0            # restore, so nothing downstream inherits the permutation

    # ---- arm B + AMBER
    if not do_amber:
        return rec
    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    E.wait_for_memory(1.5, f"phys-{pdb}")
    t0 = time.time()
    try:
        r = am.refine_coords(seq, rep, avg, k_restraint=float(k), steps=0,
                             tolerance=1.0, threads=1, memo=False)
    except NotImplementedError as exc:
        rec[f"amber_k{k}"] = dict(error=str(exc)[:120])
        return rec
    ca = np.asarray(r["ca"], float)
    bb = {a: np.asarray(r["backbone"][a], float) for a in ATOMS if a in r["backbone"]}
    phi2, psi2 = torsions_of(bb)
    rec[f"amber_k{k}"] = dict(
        rmsd=float(I.ca_rmsd(ca, nat)), energy=float(r["energy"]),
        energy_initial=float(r["energy_initial"]),
        restraint_rmsd=float(r["restraint_rmsd"]),
        moved_ca=float(I.ca_rmsd(ca, avg["CA"])),
        rama_ok=rama_ok(phi2, psi2), geom=geom_of(bb, seq),
        wall=round(time.time() - t0, 2))
    return rec


# ------------------------------------------------------------------ drivers
def _path(mode, draw):
    return os.path.join(RESULTS, f"phys_repl_{mode}{int(draw)}.json")


def run_draw(mode, draw, targets=None, k=K, resume=True, do_proj=True, do_amber=True):
    tg = targets if targets is not None else I.targets()
    p = _path(mode, draw)
    rows, done = [], set()
    if resume and os.path.exists(p):
        with open(p) as fh:
            rows = json.load(fh).get("per_target", [])
        done = {r["pdb"] for r in rows}
        if len(done) >= len(tg):
            print(f"[{mode}{draw}] already complete ({len(done)} targets)", flush=True)
            return rows
    t0 = time.time()
    for i, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        rows.append(run_target(t, mode, draw, k, do_proj=do_proj, do_amber=do_amber))
        r = rows[-1]
        a = r.get(f"amber_k{k}", {})
        pa = r.get("projected_rmsd", float("nan"))
        print(f"[{mode}{draw} {len(rows)}/{len(tg)}] {r['pdb']} proj {pa:.3f} "
              f"amber {a.get('rmsd', float('nan')):.3f} "
              f"d {a.get('rmsd', float('nan')) - pa:+.4f} "
              f"({time.time()-t0:.0f}s)", flush=True)
        if len(rows) % 10 == 0 or len(rows) == len(tg):
            _write(mode, draw, rows, len(tg), k)
    _write(mode, draw, rows, len(tg), k)
    return rows


def _write(mode, draw, rows, n_expected, k):
    with open(_path(mode, draw), "w") as fh:
        json.dump(dict(what="PHYS: replication of the k=30 AMBER relaxation effect",
                       mode=mode, draw=int(draw), k=k, per_target=rows,
                       n_rows=len(rows), n_expected=int(n_expected),
                       complete=len(rows) >= int(n_expected)), fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))


# ------------------------------------------------------------------ determinism audit
def determinism(n=6, k=K):
    """Is the k = 30 pipeline reproducible ACROSS PROCESSES?  Writes a digest per target.

    Run twice from two separate interpreters and compare.  This is the check A2 says the
    sprint never had: every earlier determinism check ran inside one process.
    """
    import hashlib
    tg = I.targets()[:n]
    out = {}
    for t in tg:
        r = run_target(t, "canon", 0, k)
        a = r[f"amber_k{k}"]
        blob = f"{r['projected_rmsd']!r}|{r['raw_avg_rmsd']!r}|{a['rmsd']!r}|{a['energy']!r}"
        out[r["pdb"]] = dict(digest=hashlib.blake2b(blob.encode(), digest_size=8).hexdigest(),
                             projected=r["projected_rmsd"], amber=a["rmsd"], energy=a["energy"])
        print(f"{r['pdb']} {out[r['pdb']]['digest']} proj {r['projected_rmsd']!r} "
              f"amber {a['rmsd']!r} E {a['energy']!r}", flush=True)
    with open(os.path.join(RESULTS, "phys_determinism_%d.json" % os.getpid()), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


# ------------------------------------------------------------------ analysis
def _load(mode, draw):
    p = _path(mode, draw)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        d = json.load(fh)
    return d["per_target"]


def _arms(rows, k=K):
    key = f"amber_k{k}"
    rows = [r for r in rows if key in r and "rmsd" in r[key]]
    rows.sort(key=lambda r: r["pdb"])
    return rows


def effect(rows, k=K):
    rows = _arms(rows, k)
    A = np.array([r["projected_rmsd"] for r in rows])
    V = np.array([r[f"amber_k{k}"]["rmsd"] for r in rows])
    folds = [r["fold"] for r in rows]
    p = I.paired(V, A, folds=folds)
    p["mean_A"] = float(A.mean()); p["mean_V"] = float(V.mean())
    p["pdbs"] = [r["pdb"] for r in rows]
    p["diff"] = (V - A).tolist()
    return p


def conc_verdict(d):
    """One PASS/FAIL verdict from the null-calibrated concentration check."""
    d = np.asarray(d, float)
    nl = concentration_null(d, seed=0)
    order = np.argsort(d)
    dt20 = float(d[order[20:]].mean()) if len(d) > 20 else float("nan")
    share = float(d[order[:10]].sum() / d.sum()) if d.sum() != 0 else float("nan")
    p20 = float((nl["dt20_draws"] < dt20).mean())
    ps = float((nl["share_draws"] < share).mean())
    conc = (p20 > 0.90) or (ps > 0.90)
    return dict(dt20=dt20, share=share, pct_dt20=p20, pct_share=ps,
                verdict="FAIL(concentrated)" if conc else "PASS",
                mean_sd=float(d.mean() / d.std(ddof=1)))


def _align(rowsA, rowsB, k=K):
    """Two draws on the common target set, in a common order."""
    a = {r["pdb"]: r for r in _arms(rowsA, k)}
    b = {r["pdb"]: r for r in _arms(rowsB, k)}
    com = sorted(set(a) & set(b))
    return [a[p] for p in com], [b[p] for p in com], com


def nullstruct(k=K, draw=1):
    """What the EXACT NULL actually found: a heavy tail, and a convergence defect.

    The frame null is zero by construction, so every non-zero entry is a defect of the
    machinery.  Two things fall out and they point in OPPOSITE directions, so both are
    printed together rather than either being quotable alone:

      1. `refine_coords` does not always converge.  Four targets end a minimisation above
         1000 kcal/mol -- 1MF6 reaches 8.9e8 kcal/mol in one frame, with two clashes and a
         0.13 geometry deviation -- and their RMSDs were silently scored into the published
         mean.  There is no convergence gate anywhere in the pipeline.  Those four alone
         move the null's mean from -0.0005 to +0.0117 and its per-target sd from 0.024 to
         0.137.
      2. Gate them out and the null collapses, but the effect turns out to live entirely on
         the targets where the minimiser is frame-SENSITIVE.  On the frame-reproducible
         majority it is not distinguishable from zero and loses more often than it wins.
    """
    can = _arms(_load("canon", 0) or [], k)
    fr = _load("frame", draw)
    if fr is None:
        print("no frame draw"); return
    fr = {r["pdb"]: r for r in _arms(fr, k)}
    com = [r for r in can if r["pdb"] in fr]
    key = f"amber_k{k}"
    V0 = np.array([r[key]["rmsd"] for r in com]); A0 = np.array([r["projected_rmsd"] for r in com])
    V1 = np.array([fr[r["pdb"]][key]["rmsd"] for r in com])
    A1 = np.array([fr[r["pdb"]].get("projected_rmsd", r["projected_rmsd"]) for r in com])
    E0 = np.array([r[key]["energy"] for r in com]); E1 = np.array([fr[r["pdb"]][key]["energy"] for r in com])
    pdbs = [r["pdb"] for r in com]
    bad = (E0 > 1e3) | (E1 > 1e3)
    print(f"=== EXACT NULL STRUCTURE (frame {draw} vs frame 0), n={len(com)} ===")
    print(f"non-converged minimisations (final E > 1e3 kcal/mol in either frame): "
          f"{[pdbs[i] for i in np.where(bad)[0]]}")
    print(f"  their final energies: canon "
          f"{['%.4g' % E0[i] for i in np.where(bad)[0]]}  frame{draw} "
          f"{['%.4g' % E1[i] for i in np.where(bad)[0]]}")
    d = V1 - V0
    print(f"\nnull |d| exceedance: " + "  ".join(
        f">{t:g}:{int((np.abs(d) > t).sum())}" for t in (1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 1.0)))
    o = np.argsort(-np.abs(d))
    print(f"largest: " + ", ".join(f"{pdbs[i]} {d[i]:+.4f}" for i in o[:6]))

    def row(nm, m):
        p0 = I.paired(V0[m], A0[m]); p1 = I.paired(V1[m], A1[m]); dn = I.paired(V1[m], V0[m])
        print(f"{nm:34s} n={int(m.sum()):3d}  effect f0 {p0['mean_diff']:+.5f} "
              f"[{p0['ci95'][0]:+.5f},{p0['ci95'][1]:+.5f}] {p0['n_better']}W/{p0['n_worse']}L"
              f" | f{draw} {p1['mean_diff']:+.5f} [{p1['ci95'][0]:+.5f},{p1['ci95'][1]:+.5f}]"
              f" | NULL {dn['mean_diff']:+.6f} [{dn['ci95'][0]:+.6f},{dn['ci95'][1]:+.6f}]"
              f" sd {d[m].std(ddof=1):.5f}")
        return dict(name=nm, n=int(m.sum()), f0=p0, f1=p1, null=dn)

    res = [row("ALL (as published)", np.ones(len(com), bool)),
           row("converged only", ~bad),
           row("  converged, minimiser STABLE", (~bad) & (np.abs(E1 - E0) <= 0.5)),
           row("  converged, minimiser UNSTABLE", (~bad) & (np.abs(E1 - E0) > 0.5))]
    with open(os.path.join(RESULTS, f"phys_nullstruct_frame{draw}.json"), "w") as fh:
        json.dump(dict(non_converged=[pdbs[i] for i in np.where(bad)[0]], rows=res), fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return res


def _armA_rows(rows):
    rows = [r for r in rows if "projected_rmsd" in r]
    rows.sort(key=lambda r: r["pdb"])
    return rows


def report(k=K, boot_draws=(1, 2, 3, 4, 5), frame_draws=(1, 2, 3, 4, 5),
           sp_draws=(1, 2, 3)):
    """The whole deliverable, printed.  Every table row carries n, CI, W/L and a verdict."""
    key = f"amber_k{k}"
    canon = _load("canon", 0)
    if canon is None:
        print("no canon draw yet"); return
    can = _arms(canon, k)
    canA = {r["pdb"]: r["projected_rmsd"] for r in can}
    canV = {r["pdb"]: r[key]["rmsd"] for r in can}
    out = {}

    print("=" * 100)
    print("PHYS -- REPLICATION OF THE k=30 AMBER RELAXATION EFFECT")
    print("=" * 100)

    # ---------------------------------------------------------------- 0. reproduction
    s14 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "s14", "results", "ener_avgrefine.json")
    if os.path.exists(s14):
        with open(s14) as fh:
            old = {r["pdb"]: r for r in json.load(fh)["per_target"] if "rmsd" in r.get(key, {})}
        com = [r for r in can if r["pdb"] in old]
        dA = max(abs(r["projected_rmsd"] - old[r["pdb"]]["projected_rmsd"]) for r in com)
        dV = max(abs(r[key]["rmsd"] - old[r["pdb"]][key]["rmsd"]) for r in com)
        print(f"\n[0] CROSS-PROCESS REPRODUCTION of the Sprint 14 artefact, n={len(com)}")
        print(f"    max |delta| arm A (projection)   {dA:.3e} A")
        print(f"    max |delta| arm B + AMBER k=30   {dV:.3e} A")
        print(f"    VERDICT: {'BIT-IDENTICAL' if max(dA, dV) == 0.0 else 'NOT bit-identical'}"
              f" -- the pipeline contains no RNG, so there is no start draw to vary.")
        out["reproduction"] = dict(n=len(com), max_dA=dA, max_dV=dV)

    # ---------------------------------------------------------------- 1. the effect, per draw
    print(f"\n[1] THE EFFECT PER DRAW  (arm B + AMBER k=30) minus (arm A projection); "
          f"negative = AMBER better")
    print(f"{'draw':22s} {'n':>4s} {'armA':>7s} {'armB+A':>7s} {'mean':>9s} {'CI95':>20s} "
          f"{'median':>8s} {'W/L':>8s} {'mean/sd':>8s} {'conc':>17s}")
    rows_tab = []

    def _line(name, V, A, folds, pdbs):
        p = I.paired(np.asarray(V), np.asarray(A), folds=folds)
        c = conc_verdict(np.asarray(V) - np.asarray(A))
        print(f"{name:22s} {p['n']:4d} {np.mean(A):7.3f} {np.mean(V):7.3f} "
              f"{p['mean_diff']:+9.4f} [{p['ci95'][0]:+.4f},{p['ci95'][1]:+.4f}] "
              f"{p['median_diff']:+8.4f} {p['n_better']:3d}/{p['n_worse']:<3d} "
              f"{c['mean_sd']:+8.3f} {c['verdict']:>17s}")
        rows_tab.append(dict(name=name, paired=p, conc=c,
                             per_fold=p.get("per_fold")))
        return p

    p0 = _line("canon (draw 0)", [r[key]["rmsd"] for r in can],
               [r["projected_rmsd"] for r in can], [r["fold"] for r in can],
               [r["pdb"] for r in can])
    boot_eff = []
    for d in boot_draws:
        rw = _load("boot", d)
        if rw is None:
            continue
        rw = _arms(rw, k)
        p = _line(f"boot draw {d}", [r[key]["rmsd"] for r in rw],
                  [r["projected_rmsd"] for r in rw], [r["fold"] for r in rw],
                  [r["pdb"] for r in rw])
        boot_eff.append(p["mean_diff"])
    frame_eff = []
    for d in frame_draws:
        rw = _load("frame", d)
        if rw is None:
            continue
        rw = _arms(rw, k)
        A = [r.get("projected_rmsd", canA[r["pdb"]]) for r in rw]
        p = _line(f"frame draw {d}", [r[key]["rmsd"] for r in rw], A,
                  [r["fold"] for r in rw], [r["pdb"] for r in rw])
        frame_eff.append(p["mean_diff"])
    sp_eff = []
    for d in sp_draws:
        rw = _load("startperm", d)
        if rw is None:
            continue
        rw = _armA_rows(rw)
        p = _line(f"startperm draw {d}", [canV[r["pdb"]] for r in rw],
                  [r["projected_rmsd"] for r in rw], [r["fold"] for r in rw],
                  [r["pdb"] for r in rw])
        sp_eff.append(p["mean_diff"])
    if boot_eff:
        print(f"\n    boot draws:  mean {np.mean(boot_eff):+.4f}  sd {np.std(boot_eff, ddof=1) if len(boot_eff) > 1 else float('nan'):.4f}"
              f"  range [{min(boot_eff):+.4f}, {max(boot_eff):+.4f}]  n_draws={len(boot_eff)}")
    if frame_eff:
        print(f"    frame draws: mean {np.mean(frame_eff):+.4f}  sd {np.std(frame_eff, ddof=1) if len(frame_eff) > 1 else float('nan'):.4f}"
              f"  range [{min(frame_eff):+.4f}, {max(frame_eff):+.4f}]  n_draws={len(frame_eff)}")
    if sp_eff:
        print(f"    startperm:   mean {np.mean(sp_eff):+.4f}  sd {np.std(sp_eff, ddof=1) if len(sp_eff) > 1 else float('nan'):.4f}"
              f"  range [{min(sp_eff):+.4f}, {max(sp_eff):+.4f}]  n_draws={len(sp_eff)}")
    out["effect_per_draw"] = rows_tab
    out["boot_effects"] = boot_eff
    out["frame_effects"] = frame_eff
    out["startperm_effects"] = sp_eff

    # ---------------------------------------------------------------- 2. the EXACT null
    print(f"\n[2] THE EXACT NULL -- the same structure relaxed in a different LAB FRAME.")
    print(f"    ff14SB/GBn2, the positional restraint and every RMSD are rigid-invariant, so")
    print(f"    each row below is ZERO BY CONSTRUCTION.  Anything it shows is numerical floor.")
    print(f"{'null pair':22s} {'n':>4s} {'mean':>9s} {'CI95':>20s} {'median':>9s} "
          f"{'W/L':>8s} {'sd':>8s} {'max|d|':>8s} {'conc':>17s}")
    nulls = []
    for d in frame_draws:
        rw = _load("frame", d)
        if rw is None:
            continue
        rw = _arms(rw, k)
        V = np.array([r[key]["rmsd"] for r in rw])
        V0 = np.array([canV[r["pdb"]] for r in rw])
        p = I.paired(V, V0, folds=[r["fold"] for r in rw])
        c = conc_verdict(V - V0)
        print(f"{'AMBER frame ' + str(d) + ' vs 0':22s} {p['n']:4d} {p['mean_diff']:+9.5f} "
              f"[{p['ci95'][0]:+.5f},{p['ci95'][1]:+.5f}] {p['median_diff']:+9.5f} "
              f"{p['n_better']:3d}/{p['n_worse']:<3d} {np.std(V - V0, ddof=1):8.5f} "
              f"{np.abs(V - V0).max():8.4f} {c['verdict']:>17s}")
        nulls.append(dict(draw=d, paired=p, conc=c, sd=float(np.std(V - V0, ddof=1)),
                          maxabs=float(np.abs(V - V0).max())))
        # arm A's own frame sensitivity, where it was computed
        hasA = [r for r in rw if "projected_rmsd" in r]
        if hasA:
            dd = np.array([r["projected_rmsd"] - canA[r["pdb"]] for r in hasA])
            print(f"{'  (arm A frame ' + str(d) + ')':22s} {len(dd):4d} {dd.mean():+9.5f} "
                  f"{'':20s} {np.median(dd):+9.5f} {'':8s} {np.std(dd, ddof=1):8.5f} "
                  f"{np.abs(dd).max():8.4f}")
            nulls[-1]["armA_sd"] = float(np.std(dd, ddof=1))
            nulls[-1]["armA_maxabs"] = float(np.abs(dd).max())
    for d in sp_draws:
        rw = _load("startperm", d)
        if rw is None:
            continue
        rw = _armA_rows(rw)
        A1 = np.array([r["projected_rmsd"] for r in rw])
        A0 = np.array([canA[r["pdb"]] for r in rw])
        p = I.paired(A1, A0, folds=[r["fold"] for r in rw])
        c = conc_verdict(A1 - A0)
        print(f"{'armA startperm ' + str(d):22s} {p['n']:4d} {p['mean_diff']:+9.5f} "
              f"[{p['ci95'][0]:+.5f},{p['ci95'][1]:+.5f}] {p['median_diff']:+9.5f} "
              f"{p['n_better']:3d}/{p['n_worse']:<3d} {np.std(A1 - A0, ddof=1):8.5f} "
              f"{np.abs(A1 - A0).max():8.4f} {c['verdict']:>17s}")
        print(f"{'':24s} basin flips (|d| > 1e-6): "
              f"{int((np.abs(A1 - A0) > 1e-6).sum())}/{len(A1)} targets")
        nulls.append(dict(kind="startperm", draw=d, paired=p, conc=c,
                          sd=float(np.std(A1 - A0, ddof=1)),
                          maxabs=float(np.abs(A1 - A0).max()),
                          n_flip=int((np.abs(A1 - A0) > 1e-6).sum())))
    if nulls:
        fl = max(abs(x["paired"]["ci95"][0]) for x in nulls) if nulls else float("nan")
        fl = max(fl, max(abs(x["paired"]["ci95"][1]) for x in nulls))
        print(f"\n    FLOOR FOR THIS PIPELINE: widest exact-null CI bound = {fl:.5f} A"
              f"   (the RESTRAINT workstream's floor, from a different pipeline, was 0.081 A)")
        out["floor"] = float(fl)
    out["nulls"] = nulls

    # ------------------------------------------------- 1b. the confound the boot axis carries
    print(f"\n[1b] THE CONFOUND IN THE BOOT AXIS, stated because it biases the answer.")
    print(f"     A bootstrap of 75 windows keeps only ~47 distinct ones, so the average is")
    print(f"     noisier, the backbone contracts MORE, and the projection cost (A - B) rises.")
    print(f"     The effect is known to scale with that cost (rho = -0.51 on canon), so a")
    print(f"     LARGER effect on boot draws is mechanism, not evidence of stability.")
    print(f"{'draw':22s} {'n_distinct':>10s} {'proj cost A-B':>14s} {'effect':>10s}")
    ctx = []
    for nm, rw in ([("canon (draw 0)", can)] +
                   [(f"boot draw {d}", _arms(_load('boot', d) or [], k)) for d in boot_draws
                    if _load('boot', d)]):
        if not rw:
            continue
        cost = float(np.mean([r["projected_rmsd"] - r["raw_avg_rmsd"] for r in rw]))
        eff = float(np.mean([r[key]["rmsd"] - r["projected_rmsd"] for r in rw]))
        nd = float(np.mean([r.get("n_distinct", r["n_windows"]) for r in rw]))
        print(f"{nm:22s} {nd:10.1f} {cost:+14.4f} {eff:+10.4f}")
        ctx.append(dict(name=nm, n_distinct=nd, proj_cost=cost, effect=eff))
    if len(ctx) > 2:
        from scipy import stats as _st
        rr = _st.pearsonr([c["proj_cost"] for c in ctx], [c["effect"] for c in ctx])
        print(f"     across draws: r(proj cost, effect) = {rr[0]:+.3f} (n_draws={len(ctx)})")
    out["boot_context"] = ctx

    # ------------------------------------------------- 2b. per-fold and FAIL18, every draw
    print(f"\n[2b] PER-FOLD and FAIL18 detail.  (FAIL18 is reported for cross-sprint")
    print(f"     comparability only; per the Phase 0 audit it is a BAND=1.5 threshold")
    print(f"     artefact and no argument may rest on its membership.)")
    f18 = set(I.FAIL18)
    print(f"{'draw':22s} " + " ".join(f"{'f'+str(i):>8s}" for i in range(5)) +
          f" {'FAIL18':>9s} {'other':>9s}")
    for nm, rw in ([("canon (draw 0)", can)] +
                   [(f"boot draw {d}", _arms(_load('boot', d) or [], k)) for d in boot_draws
                    if _load('boot', d)]):
        if not rw:
            continue
        dd = np.array([r[key]["rmsd"] - r["projected_rmsd"] for r in rw])
        fo = np.array([r["fold"] for r in rw])
        m18 = np.array([r["pdb"] in f18 for r in rw])
        print(f"{nm:22s} " + " ".join(f"{dd[fo==i].mean():+8.4f}" if (fo == i).any() else f"{'-':>8s}"
                                      for i in range(5)) +
              f" {dd[m18].mean():+9.4f} {dd[~m18].mean():+9.4f}")

    # ---------------------------------------------------------------- 3. accuracy vs validity
    print(f"\n[3] ACCURACY versus VALIDITY -- the two things relaxation does, separated.")
    folds = [r["fold"] for r in can]

    def _val(name, v, a, fmt="{:+.4f}"):
        p = I.paired(np.asarray(v, float), np.asarray(a, float), folds=folds)
        print(f"  {name:44s} {np.mean(a):7.4f} -> {np.mean(v):7.4f}   "
              f"{fmt.format(p['mean_diff'])} [{p['ci95'][0]:+.4f},{p['ci95'][1]:+.4f}]  "
              f"{p['n_better']:3d}W/{p['n_worse']:<3d}L  median {p['median_diff']:+.4f}")
        return p

    print(f"  {'quantity':44s} {'armA':>7s}    {'armB+A':>7s}   {'paired diff (neg = AMBER lower)'}")
    out["accuracy"] = _val("ACCURACY: CA-RMSD to native (A)  [ORACLE]",
                           [r[key]["rmsd"] for r in can], [r["projected_rmsd"] for r in can])
    out["rama_matched"] = _val("VALIDITY: 1 - Ramachandran-favoured fraction",
                               [1 - r[key]["rama_ok"] for r in can],
                               [1 - r["projected_rama_ok"] for r in can])
    if "projected_lam_rama_ok" in can[0]:
        out["rama_s14"] = _val("  (the same, vs the lam=0.3 arm S14 used by mistake)",
                               [1 - r[key]["rama_ok"] for r in can],
                               [1 - r["projected_lam_rama_ok"] for r in can])
    if "projected_geom" in can[0]:
        out["clash"] = _val("VALIDITY: n heavy-atom clashes < 2.0 A",
                            [r[key]["geom"]["n_clash_2A"] for r in can],
                            [r["projected_geom"]["n_clash_2A"] for r in can])
        out["minheavy"] = _val("VALIDITY: min heavy-atom separation (A, higher better)",
                               [r[key]["geom"]["min_heavy"] for r in can],
                               [r["projected_geom"]["min_heavy"] for r in can])
    g = np.array([r[key]["geom"]["geom_rms_rel_dev"] for r in can])
    print(f"  {'GEOMETRY: rms relative bond/angle deviation':44s} "
          f"{0.0:7.4f} -> {g.mean():7.4f}   (arm A is ideal by construction; max {g.max():.4f})")
    out["geom_dev"] = dict(mean=float(g.mean()), max=float(g.max()))
    print(f"  {'RAMA fractions (higher better):':44s} "
          f"armA lam=0 {np.mean([r['projected_rama_ok'] for r in can]):.4f}   "
          f"armA lam=0.3 {np.mean([r.get('projected_lam_rama_ok', np.nan) for r in can]):.4f}   "
          f"raw average {np.mean([r['raw_avg_rama_ok'] for r in can]):.4f}   "
          f"AMBER {np.mean([r[key]['rama_ok'] for r in can]):.4f}")

    # ---------------------------------------------------------------- 4. cost
    print(f"\n[4] COST ACCOUNTING")
    wa = np.array([r[key]["wall"] for r in can])
    wp = np.array([r["projected_wall"] for r in can if "projected_wall" in r])
    print(f"    arm B + AMBER k=30 (refine_coords, steps=0, tol 1.0): "
          f"mean {wa.mean():.2f} s  median {np.median(wa):.2f} s  max {wa.max():.2f} s")
    if len(wp):
        print(f"    arm A projection (lam_path, multi=True, 9 L-BFGS runs): "
              f"mean {wp.mean():.2f} s  median {np.median(wp):.2f} s  max {wp.max():.2f} s")
    print(f"    NOTE: the audited 8.3-23.3 ms figure is an AMBER SINGLE POINT.  This arm is a")
    print(f"    full restrained MINIMISATION, {wa.mean()/0.0083:.0f}x that.  No budget-matched")
    print(f"    comparison enters the -0.022 A figure, so the 28 ms -> 8-14 ms correction")
    print(f"    does not move it.")
    out["cost"] = dict(amber_wall_mean=float(wa.mean()),
                       proj_wall_mean=float(wp.mean()) if len(wp) else None)

    with open(os.path.join(RESULTS, "phys_repl_report.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return out


if __name__ == "__main__":
    a = sys.argv[1:]

    def _opt(name, default=None):
        return a[a.index(name) + 1] if name in a else default

    mode = _opt("--mode")
    if "--nullstruct" in a:
        nullstruct(draw=int(_opt("--draws", "1")))
    elif "--report" in a:
        report()
    elif mode == "determinism":
        determinism(int(_opt("--n", "6")))
    elif mode:
        draws = [int(x) for x in _opt("--draws", "0").split(",")]
        tg = I.targets()
        lim = _opt("--limit")
        if lim:
            tg = tg[:int(lim)]
        noproj = "--noproj" in a
        noamber = "--noamber" in a
        for d in draws:
            run_draw(mode, d, targets=tg, do_proj=not noproj, do_amber=not noamber)
    else:
        print(__doc__)
