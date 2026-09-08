"""s16/repair.py -- AMBER AS A REPAIR OPERATOR (the REPAIR workstream, Sprint 16).

NOBODY HAS TESTED THIS.  The ENERGY workstream tests AMBER as a RANKING FILTER over a fixed
retrieval ensemble; the coordinator tests it as a filter over VQE ensembles.  Both consume
AMBER as a scalar.  The user's mandated architecture ends with "AMBER all-atom
repair/refinement" -- an OPERATOR THAT CHANGES COORDINATES -- and that operator has never been
swept, never been frame-nulled at any setting but one, and its DISPLACEMENT has never been
looked at.

TWO QUESTIONS.

Q1  Does restrained AMBER minimisation repair, and at what cost?  Sweep the restraint constant
    k and the minimisation depth, n = 126, convergence gate applied and excluded count printed
    every time, ACCURACY and VALIDITY as separate axes, and the rotated-lab-frame null -- which
    is ZERO BY RIGID INVARIANCE -- re-measured at EVERY setting.  An arm whose frame null is
    non-zero is not reportable.

Q2  Repair versus the programme's new geometric fact.  `s16/STEER_FINDINGS.md` establishes that
    RMSD is strongly CONVEX in torsion displacement toward the truth (half-way along the EXACT
    native error buys 7.7% [-0.0, +14.9] and makes 49/126 targets worse), and `s16/csteer.py`
    that from the pool coordinate average every native-free direction has cosine with the true
    residual indistinguishable from zero.  So: AMBER MOVES THE STRUCTURE -- IN WHICH DIRECTION,
    AND HOW FAR?  Measure the cosine of its displacement with the true residual (ORACLE), the
    magnitude, and whether it lands in the loud or the quiet half of the superposed-CA Jacobian
    spectrum -- and the same three for the ideal-geometry projection, so the two repair
    operators are comparable.

WHAT IS ALREADY KNOWN AND IS NOT RE-LITIGATED HERE (`s15/phys_FINDINGS.md`):
  * the k = 30 accuracy effect is -0.022 [-0.036, -0.009] ungated, and is ABSENT on the 84
    targets (69%) where the minimiser is frame-reproducible: -0.009 [-0.024, +0.006], 38W/46L,
    replicated on three frames;
  * the ungated frame null is +0.0117 and collapses to -0.0005 [-0.0046, +0.0035] gated;
  * `refine_coords` had no convergence gate (fixed this sprint in `core.amber`);
  * the VALIDITY half is large: Rama-favoured 0.466 -> 0.874 (116W/2L), clashes 1.397 -> 0.000
    (63W/0L).
This module asks whether ANY (k, depth) rescues the ACCURACY half, and what the operator is
doing to the geometry.  It does not re-run the k = 30 canon; it reads arm A from
`s15/results/phys_repl_canon0.json` so that every comparison is against the SAME projection
PHYS measured, and re-derives arm A's CA coordinates for Q2.

CONTROLS CARRIED IN EVERY CLAIM (AUDIT's finding: this programme's evidence read ~5x stronger
than it was because arms lacked exactly these two).

  ZERO-INFORMATION references
    `none`   do nothing -- the raw all-atom coordinate average.  A repair operator that uses no
             information at all.  It is the BEST of the three on accuracy (3.050 A) and the
             worst on validity, which is the whole point.
    `helix`  a constant ideal alpha-helix (phi = -57, psi = -47) of the same length, superposed
             on the input.  ZERO information about the target beyond its length.  It has
             perfect Ramachandran statistics and no clashes BY CONSTRUCTION, so it prices the
             VALIDITY axis: any validity claim a constant helix reproduces is not a result.
  RANDOM-DIRECTION control at MATCHED MAGNITUDE
    `rand`   the input CA displaced by an isotropic random direction scaled to the SAME RMS
             magnitude the operator moved, 3 draws under `stable_rng`.  It prices the ACCURACY
             axis: a move of that size in no particular direction.

SEEDING  `s15.seed.stable_rng` everywhere; no `hash()`.  The AMBER path has no RNG and
reproduces bit-identically across interpreters (PHYS S1), so the 0.08 A multi-start
false-positive floor DOES NOT transfer to it; its own floor is ~0.004 A gated and up to
0.038 A ungated.

COST  Every arm here is a FULL RESTRAINED MINIMISATION, ~9-19 s each on this box at
threads = 1 (measured, this module reports its own wall clock), i.e. ~1,500x an AMBER single
point and ~2-3x the ideal-geometry projection.  A single-point cost is NEVER quoted for these
arms.

    python -m s16.repair --mode sweep            # Q1, n = 126, both frames  (~3.5 h)
    python -m s16.repair --mode verify           # ladder-vs-refine_coords bit check
    python -m s16.repair --mode report           # every table
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import argparse
import json
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                        # noqa: E402
from s15.seed import stable_rng                        # noqa: E402
from s14.avgspace import top75_windows                 # noqa: E402
from s14.ener_avgrefine import ATOMS, geom_of, rama_ok, torsions_of   # noqa: E402
from s15.phys_repl import averaged_backbone_from, project_both        # noqa: E402
from core.amber import CONVERGE_MAX_KCAL               # noqa: E402

S15R = os.path.join(ROOT, "s15", "results")

#: THE SWEEP, declared here before any of it ran.
#: k = 0 IS the unrestrained arm -- the restraint force constant is a context parameter and
#: setting it to zero removes the restraint exactly.  It is listed twice in the brief
#: ("k = 0 ... and unrestrained") and it is one arm, said so plainly.
K_LADDER = (0.0, 5.0, 15.0, 30.0, 60.0)
#: minimisation depth = OpenMM `maxIterations`; 0 means "until converged", the incumbent.
DEPTH_LADDER = (25, 100, 400)
K_FOR_DEPTH = 30.0
#: frame-reproducibility band, pre-declared, native-free: |E(frame B) - E(frame A)| <= 0.5.
#: This is PHYS's definition, reused unchanged so the strata are comparable across sprints.
STABLE_DE = 0.5
N_RAND = 3
#: RADIANS.  `core.geometry.build_backbone` and `s14.ener_avgrefine.torsions_of` both work in
#: radians; the first draft of this module wrote -57.0 (degrees) here and built a nonsense
#: helix.  Caught before any number was reported; the affected rows in `repair_sweep.json`
#: are recomputed by `--mode fixhelix`, which touches no AMBER arm.
HELIX_PHI, HELIX_PSI = float(np.radians(-57.0)), float(np.radians(-47.0))


def _settings():
    out = [("k%g_full" % k, float(k), 0) for k in K_LADDER]
    out += [("k%g_d%d" % (K_FOR_DEPTH, d), float(K_FOR_DEPTH), int(d)) for d in DEPTH_LADDER]
    return out


SETTINGS = _settings()
BY_NAME = {nm: (nm, k, s) for nm, k, s in SETTINGS}

#: THE PASSES, ordered so the DECISIVE answers land first.  Re-scoped by the coordinator
#: after the first 5 targets of a single 8-setting pass showed a ~4.2 h path that was
#: starving sibling agents; the 5-target partial is preserved as
#: `results/repair_ladder5_superseded.json` and is quoted nowhere.
#:
#:   A  the shipped setting k = 30, both lab frames, plus everything Q2 needs.  This alone
#:      settles the CONVERGENCE GATE and the ROTATED-FRAME NULL at full n, and if those do
#:      not come out clean no accuracy number from any k is reportable.
#:   B  the unrestrained arm k = 0.  With A this is the contrast the record turns on.
#:   C  minimisation DEPTH at k = 30 (25 / 100 / 400 iterations against "until converged").
#:   D  the SHAPE check in k, on a PRE-REGISTERED native-free stratified subsample (see
#:      `shape_subsample`), never at full n and never used to select anything.
PASSES = {"A": ["k30_full"], "B": ["k0_full"],
          "C": ["k30_d25", "k30_d100", "k30_d400"],
          "D": ["k5_full", "k15_full", "k60_full"]}
SHAPE_N = 30


# ------------------------------------------------------------------ helpers
def _sup(X, T):
    """Superpose X onto T (proper rotations only, via the frozen instrument)."""
    return I.superpose_batch(np.asarray(X, float)[None], np.asarray(T, float))[0]


def random_rigid(rng):
    """Proper rotation (det = +1) and translation.  A reflection is NOT a symmetry of ff14SB."""
    Q, R = np.linalg.qr(rng.normal(size=(3, 3)))
    Q = Q * np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1.0
    return Q, rng.normal(scale=10.0, size=3)


def helix_backbone(n):
    """ZERO-INFORMATION reference: the constant ideal alpha-helix of length n."""
    from core import geometry as geo
    bb = geo.build_backbone(np.full(n, HELIX_PHI), np.full(n, HELIX_PSI))
    return {a: np.asarray(bb[a], float) for a in ATOMS if a in bb}


def _cached_armA():
    """Arm A (the ideal-geometry projection) exactly as PHYS measured it, frame 0."""
    p = os.path.join(S15R, "phys_repl_canon0.json")
    if not os.path.exists(p):
        return {}
    with open(p) as fh:
        rows = json.load(fh)["per_target"]
    return {r["pdb"]: r for r in rows}


# ------------------------------------------------------------------ one target
def shape_subsample(n=SHAPE_N):
    """PRE-REGISTERED, NATIVE-FREE subsample for the pass-D shape check in k.

    Stratified by FOLD (5) and by CHAIN LENGTH (below / above the median length), drawn
    under `stable_rng` from the target list alone.  It reads no RMSD, no energy and no
    native coordinate, and it is fixed by this function -- printing it before pass D runs
    is what makes it pre-registered rather than chosen.  Nothing is selected on it: pass D
    is a SHAPE check and any k it favours is reported as such, never re-quoted at full n as
    an out-of-sample choice.
    """
    tg = I.targets()
    L = np.array([len(t["seq"]) for t in tg])
    med = float(np.median(L))
    strata = {}
    for t, ln in zip(tg, L):
        strata.setdefault((int(t["fold"]), bool(ln > med)), []).append(t["pdb"])
    keys = sorted(strata)
    rng = stable_rng("s16repair", "shape_subsample", int(n))
    per = max(1, int(round(n / len(keys))))
    out = []
    for k in keys:
        pool = sorted(strata[k])
        take = min(per, len(pool))
        idx = rng.permutation(len(pool))[:take]
        out += [pool[i] for i in sorted(idx)]
    out = sorted(set(out))
    return out


def run_target(t, armA_cache, frames=(0, 1), do_proj=True, settings=None):
    """One target: the requested settings, both lab frames.  One OpenMM builder throughout."""
    from core import amber as am
    from core import geometry as geo
    import torsion_lib2 as tl2

    pdb, seq, fold = t["pdb"], t["seq"], int(t["fold"])
    n = len(seq)
    W, PHI, PSI, u = top75_windows(pdb)
    nat = np.asarray(u["nat_ca"], float)
    avg, C_ca, dev = averaged_backbone_from(W, PHI, PSI)

    rec = dict(pdb=pdb, fold=fold, n=n, ca_identity_max_dev=float(dev),
               raw_avg_rmsd=float(I.ca_rmsd(avg["CA"], nat)),
               raw_avg_rama_ok=float(rama_ok(*torsions_of(avg))),
               raw_avg_geom=geom_of(avg, seq),
               input_ca=np.asarray(avg["CA"], float).tolist(),
               nat_ca=nat.tolist())

    # ---- ZERO-INFORMATION reference: the constant alpha-helix, superposed on the input.
    hb = helix_backbone(n)
    hca = _sup(hb["CA"], avg["CA"])
    rec["helix"] = dict(rmsd=float(I.ca_rmsd(hb["CA"], nat)),
                        rama_ok=float(rama_ok(np.full(n, HELIX_PHI), np.full(n, HELIX_PSI))),
                        geom=geom_of(hb, seq),
                        moved_ca=float(I.ca_rmsd(hca, avg["CA"])),
                        ca=hca.tolist())

    # ---- arm A, the ideal-geometry projection.  Coordinates kept for Q2.
    if do_proj:
        t0 = time.time()
        pr = project_both(C_ca, seq, fold)
        bbA = geo.build_backbone(pr["fit_phi"], pr["fit_psi"])
        rec["proj"] = dict(rmsd=float(I.ca_rmsd(pr["fit_ca"], nat)),
                           rama_ok=float(rama_ok(pr["fit_phi"], pr["fit_psi"])),
                           geom=geom_of({a: np.asarray(bbA[a], float)
                                         for a in ATOMS if a in bbA}, seq),
                           moved_ca=float(I.ca_rmsd(pr["fit_ca"], avg["CA"])),
                           ca=np.asarray(pr["fit_ca"], float).tolist(),
                           phi=np.asarray(pr["fit_phi"], float).tolist(),
                           psi=np.asarray(pr["fit_psi"], float).tolist(),
                           wall=round(time.time() - t0, 2))
        c = armA_cache.get(pdb)
        if c is not None:
            rec["proj"]["cached_rmsd"] = float(c["projected_rmsd"])
            rec["proj"]["cache_delta"] = float(rec["proj"]["rmsd"] - c["projected_rmsd"])

    # ---- AMBER, every setting, every frame.
    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    for fi in frames:
        if fi == 0:
            coords = avg
        else:
            Q, tv = random_rigid(stable_rng("s16repair", "frame", int(fi), pdb))
            coords = {a: v @ Q.T + tv for a, v in avg.items()}
        for name, k, steps in (settings if settings is not None else SETTINGS):
            t0 = time.time()
            r = am.refine_coords(seq, rep, coords, k_restraint=float(k), steps=int(steps),
                                 tolerance=1.0, threads=1, memo=False)
            ca = np.asarray(r["ca"], float)
            bb = {a: np.asarray(r["backbone"][a], float) for a in ATOMS if a in r["backbone"]}
            ph2, ps2 = torsions_of(bb)
            d = dict(rmsd=float(I.ca_rmsd(ca, nat)),
                     energy=float(r["energy"]), energy_initial=float(r["energy_initial"]),
                     converged=bool(r["converged"]),
                     restraint_rmsd=float(r["restraint_rmsd"]),
                     moved_ca=float(I.ca_rmsd(ca, coords["CA"])),
                     rama_ok=float(rama_ok(ph2, ps2)),
                     geom=geom_of(bb, seq),
                     wall=round(time.time() - t0, 2))
            if fi == 0:
                d["ca"] = ca.tolist()
                d["phi"] = np.asarray(ph2, float).tolist()
                d["psi"] = np.asarray(ps2, float).tolist()
            rec[f"f{fi}_{name}"] = d
    return rec


# ------------------------------------------------------------------ driver
def _path(tag):
    return os.path.join(RESULTS, f"repair_{tag}.json")


def run_sweep(tag="A", targets=None, resume=True, frames=(0, 1), settings=None,
              do_proj=True):
    tg = targets if targets is not None else I.targets()
    settings = settings if settings is not None else SETTINGS
    p = _path(tag)
    rows, done = [], set()
    if resume and os.path.exists(p):
        with open(p) as fh:
            rows = json.load(fh).get("per_target", [])
        done = {r["pdb"] for r in rows}
    cache = _cached_armA()
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        rows.append(run_target(t, cache, frames=frames, do_proj=do_proj,
                               settings=settings))
        r = rows[-1]
        msg = " ".join(f"{nm}={r[f'f0_{nm}']['rmsd']:.3f}" for nm, _, _ in settings)
        print(f"[{len(rows)}/{len(tg)}] {r['pdb']} n={r['n']} raw={r['raw_avg_rmsd']:.3f} "
              f"proj={r.get('proj',{}).get('rmsd',float('nan')):.3f} {msg} "
              f"({time.time()-t0:.0f}s)", flush=True)
        if len(rows) % 5 == 0 or len(rows) == len(tg):
            _write(tag, rows, len(tg), settings)
    _write(tag, rows, len(tg), settings)
    return rows


def _write(tag, rows, n_expected, settings=SETTINGS):
    with open(_path(tag), "w") as fh:
        json.dump(dict(what="REPAIR: AMBER as an operator that changes coordinates",
                       settings=[list(s) for s in settings],
                       converge_max_kcal=float(CONVERGE_MAX_KCAL),
                       per_target=rows, n_rows=len(rows), n_expected=int(n_expected)),
                  fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))


# ------------------------------------------------------------------ verification
def verify(n=6):
    """The depth arms use `refine_coords(steps=N)`, the same call the incumbent uses with
    N = 0.  This checks that this module's `full` arm reproduces the PHYS artefact
    BIT-IDENTICALLY, which is the only claim that matters for comparability."""
    cache = _cached_armA()
    tg = I.targets()[:n]
    out = []
    for t in tg:
        r = run_target(t, cache, frames=(0,), do_proj=True)
        c = cache.get(t["pdb"], {})
        a = c.get("amber_k30.0", {})
        row = dict(pdb=t["pdb"],
                   mine=r["f0_k30_full"]["rmsd"], phys=a.get("rmsd"),
                   d_amber=(r["f0_k30_full"]["rmsd"] - a["rmsd"]) if a else None,
                   mine_E=r["f0_k30_full"]["energy"], phys_E=a.get("energy"),
                   d_proj=r["proj"].get("cache_delta"))
        out.append(row)
        print(f"{row['pdb']} amber d={row['d_amber']:+.3e} proj d={row['d_proj']:+.3e} "
              f"E d={row['mine_E'] - a['energy']:+.3e}", flush=True)
    with open(_path("verify"), "w") as fh:
        json.dump(out, fh, indent=1)
    mx = max(abs(r["d_amber"]) for r in out if r["d_amber"] is not None)
    mp = max(abs(r["d_proj"]) for r in out if r["d_proj"] is not None)
    print(f"max |d| amber {mx:.3e}   projection {mp:.3e}   "
          f"({'BIT-IDENTICAL' if mx == 0 and mp == 0 else 'DIVERGENT'})")
    return out


def fixhelix(tag="sweep"):
    """Recompute the ZERO-INFORMATION helix reference in place.

    The first draft built the constant helix from DEGREES where the builder wants RADIANS.
    No AMBER arm and no accuracy number depends on it, so this repairs the artefact without
    re-running a single minimisation; the old block is kept as `helix_degbug` so the error is
    preserved in place rather than edited away.
    """
    p = _path(tag)
    with open(p) as fh:
        d = json.load(fh)
    tg = {t["pdb"]: t for t in I.targets()}
    for r in d["per_target"]:
        t = tg[r["pdb"]]
        seq, n = t["seq"], r["n"]
        A = np.asarray(r["input_ca"], float)
        nat = np.asarray(r["nat_ca"], float)
        hb = helix_backbone(n)
        hca = _sup(hb["CA"], A)
        r["helix_degbug"] = r.get("helix")
        r["helix"] = dict(rmsd=float(I.ca_rmsd(hb["CA"], nat)),
                          rama_ok=float(rama_ok(np.full(n, HELIX_PHI), np.full(n, HELIX_PSI))),
                          geom=geom_of(hb, seq),
                          moved_ca=float(I.ca_rmsd(hca, A)),
                          ca=hca.tolist())
    with open(p, "w") as fh:
        json.dump(d, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    H = [r["helix"] for r in d["per_target"]]
    print("helix (ZERO-INFORMATION) recomputed on %d rows: rmsd %.3f  rama %.3f  clash %.3f"
          % (len(H), np.mean([h["rmsd"] for h in H]), np.mean([h["rama_ok"] for h in H]),
             np.mean([h["geom"]["n_clash_2A"] for h in H])))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="pass",
                    choices=("pass", "verify", "fixhelix", "prereg"))
    ap.add_argument("--pass", dest="pas", default="A", choices=sorted(PASSES))
    ap.add_argument("--tag", default=None)
    ap.add_argument("--n", type=int, default=0)
    a = ap.parse_args()
    if a.mode == "verify":
        verify(a.n or 6)
    elif a.mode == "fixhelix":
        fixhelix(a.tag or "A")
    elif a.mode == "prereg":
        sub = shape_subsample()
        print("PRE-REGISTERED shape subsample for pass D, n = %d (native-free, "
              "fold x length stratified, stable_rng):" % len(sub))
        print("  " + " ".join(sub))
    else:
        st = [BY_NAME[nm] for nm in PASSES[a.pas]]
        tg = I.targets()
        if a.pas == "D":
            sub = set(shape_subsample())
            tg = [t for t in tg if t["pdb"] in sub]
            print("pass D on the PRE-REGISTERED subsample, n = %d" % len(tg), flush=True)
        if a.n:
            tg = tg[:a.n]
        # only pass A pays for the ideal-geometry projection; later passes merge onto it.
        run_sweep(tag=a.tag or a.pas, targets=tg, settings=st, do_proj=(a.pas == "A"))
