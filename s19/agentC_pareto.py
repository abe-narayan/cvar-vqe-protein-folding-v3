"""s19/agentC_pareto.py -- Q2: AMBER AS A CONSTRAINED TERMINAL OPERATOR.  THE FRONTIER.

Pre-registration: `s19/PREREG_C.md` section 4, written before this module produced a number.

WHAT THIS IS NOT.  It is not an attempt to make AMBER a stronger minimiser, and it does not
reopen AMBER as a ranker (closed, `s19/BRIEF.md` section 4).  AMBER's standing role --
stereochemical repair of the coordinate average, at a measured accuracy cost of
+0.133 A [+0.112, +0.165] at k = 30 and an irreducible +0.164 A repair tax -- is the premise, not
the target.

THE QUESTION.  The repair is bought with DISPLACEMENT: the k = 30 restraint lets every restrained
atom move, and the Ca trace the metric scores moves with it.  Is there a differently-CONSTRAINED
AMBER that buys the same stereochemistry for less Ca displacement?  That is a two-axis question and
it is reported as two axes:

    ACCURACY   Ca-RMSD, and its delta against the arm's OWN input
    VALIDITY   the full `s16.energy_lib.panel` vector -- clashes, min heavy distance, bond and
               angle strain, Ramachandran, cis fraction, chirality -- NEVER fused into a scalar

GENUINE AMBER THROUGHOUT.  ff14SB/GBn2 through `core.amber`; no surrogate.  The only thing any arm
changes is the POSITIONAL RESTRAINT, which is not part of ff14SB: its strength `k`, which atoms it
holds, and what it holds them to.  The force field, the solvent model, the topology, the hydrogen
placement and `LocalEnergyMinimizer` are untouched in every arm.

TWO IMPLEMENTATION FACTS, RECORDED BECAUSE THEY COULD SILENTLY CORRUPT AN ARM:

  1. `core.amber.builder_for`'s cache key does NOT include the restrained atom set, and neither
     does `_memo_tag`.  A `caonly` builder and a default builder would therefore SHARE a cache
     entry and a memo entry.  This module keeps its own builder cache keyed on the restraint set
     and runs every minimisation with `memo=False`.  Verified in gate GC2.
  2. `refine_coords` rebuilds side chains from the backbone dict on every call, so a two-stage
     protocol chained through coordinates would silently re-idealise the side chains between
     stages.  The staged arms therefore chain through FULL-SYSTEM POSITIONS (`_run_ref` below),
     which is the same arithmetic as `core.amber._run` with the restraint reference made an
     explicit argument instead of being implicitly the input.

    python -m s19.agentC_pareto --smoke     # 4 targets, writes _SMOKE_*, never a result
    python -m s19.agentC_pareto             # n = 126
"""
from __future__ import annotations

import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                                   # noqa: E402
from s14.avgspace import top75_windows                            # noqa: E402
from s15.phys_repl import averaged_backbone_from, random_rigid     # noqa: E402
from s15 import seed as SD                                        # noqa: E402
from s16 import energy_lib as EL                                  # noqa: E402
from s18 import phys_lib as PL                                    # noqa: E402
from s19 import agentC_lib as CL                                  # noqa: E402

ATOMS = ("N", "CA", "C", "O", "CB")
#: COMPUTE DEVIATION FROM THE PRE-REGISTRATION, recorded before any accuracy number from the
#: dropped arms was read.  PREREG_C.md section 4.3 declared the ladder as
#: (1, 3, 10, 30, 100, 300, 1000).  `LocalEnergyMinimizer` is called with `steps = 0`, which is
#: OpenMM's "iterate until converged" and is UNBOUNDED.  At k = 1 and k = 3 the restraint is too
#: soft to hold the fold and the minimisation does not terminate in usable time: on 1A1P the
#: soft arms had not returned after 40 minutes, against ~10 s at k >= 10.  The two softest rungs
#: are therefore DROPPED, on a timing fact observed before their Ca-RMSD or validity was looked
#: at.  Capping their iterations instead was rejected because a capped arm is not the same
#: operator as the uncapped incumbent, and the incumbent's protocol is not up for change.
#: The surviving ladder still spans two decades around k = 30 and still answers "hard vs soft".
K_LADDER = (1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0)   # the FULL pre-registered ladder
K_INCUMBENT = 30.0
#: SECOND COMPUTE DEVIATION, declared here and verified in gate GC3 before any number was read.
#: `steps = 0` is OpenMM's UNBOUNDED "iterate to convergence", and on some targets it does not
#: terminate: 1D0W spent >35 minutes of a full core inside one ladder call before this module was
#: stopped (the 5-target partial is kept as `_PARTIAL_agentC_pareto_uncapped_n5.json` and is not
#: quoted anywhere).  Every arm, INCLUDING the incumbent k = 30, is therefore run with the SAME
#: finite iteration bound.  Gate GC3 checks the bound is INERT: capped k = 30 must be
#: bit-identical to uncapped k = 30.  Any arm that hits the bound is counted and reported.
#:
#: THE BOUND WAS TIGHTENED ONCE, AND THE REASON IS A MEASUREMENT, NOT A PREFERENCE.  The first
#: value, 10000, was chosen as ~12x the 708-847 FORCE EVALUATIONS a normal call needs -- but
#: `maxIterations` counts L-BFGS ITERATIONS, each of which runs a line search over several
#: evaluations, so 10000 iterations is far more than 10000 evaluations.  Measured: one capped
#: k = 10 call on 1CEK still exceeded SEVEN MINUTES, i.e. the pathological inputs need >10000
#: iterations and the bound was not bounding anything useful.  2000 is ~4-5x the iterations a
#: converging call uses and bounds one call to a few tens of seconds.  GC3 is re-run at the new
#: value; if it does not come back bit-identical for the incumbent k = 30, the bound is NOT
#: inert and Q2 does not run.
STEPS = 0
TOL = 1.0
BLEND = (0.25, 0.50, 0.75)

#: restrained atom sets.  `("N","CA","C")` is the incumbent `core.amber.RESTRAINED_BACKBONE`.
#: DEVIATION FROM THE PRE-REGISTRATION, recorded rather than repaired silently.  PREREG_C.md
#: section 4.3 declared the wider restrained set as {N, CA, C, O, CB}.  **CB does not exist for
#: glycine**, and `core.amber._index_topology` builds `_restraint_idx` by direct dict lookup, so
#: that set raises `KeyError: (i, 'CB')` on every target containing a G -- which is most of them.
#: The arm is therefore run as {N, CA, C, O}, which is the same claim ("restrain more of the
#: backbone than N/CA/C") and is constructible on all 126 targets.  It is named `bbo`, not
#: `heavybb`, so no table can read as the arm that was declared.
SETS = {"bb": ("N", "CA", "C"), "caonly": ("CA",), "bbo": ("N", "CA", "C", "O")}

_BUILDERS = {}


def builder(seq, rep, restrained):
    """A builder whose POSITIONAL RESTRAINT holds `restrained`, cached on that set.

    `core.amber.RESTRAINED_BACKBONE` is read inside `_index_topology`, so it is patched for the
    construction only and restored immediately; nothing else in `core.amber` is touched, and the
    module constant is left exactly as found.
    """
    import core.amber as am
    key = (seq, tuple(restrained))
    h = _BUILDERS.get(key)
    if h is not None:
        return h
    old = am.RESTRAINED_BACKBONE
    am.RESTRAINED_BACKBONE = tuple(restrained)
    try:
        h = am.AmberHamiltonian(seq, rep, restraint_k=am.K_MODERATE, minimization_steps=0,
                                platform_name="CPU", collapse_floor=float("-inf"), threads=1)
    finally:
        am.RESTRAINED_BACKBONE = old
    assert len(h._restraint_idx) == len(seq) * len(restrained), "restraint set did not take"
    _BUILDERS[key] = h
    return h


def drop_builders():
    import core.amber as am
    for h in list(_BUILDERS.values()):
        am._drop_builder(h)
    _BUILDERS.clear()


def _run_ref(H, pos_nm, ref_nm, k_restraint, steps=STEPS, tolerance=TOL):
    """`core.amber._run`, with the restraint REFERENCE as an explicit argument.

    Identical arithmetic and identical reporting; the single change is that the harmonic
    restraint is anchored at `ref_nm` instead of at the input positions, which is what a
    pull-back arm needs and what chaining two stages through full-system positions needs.
    `ref_nm is pos_nm` reproduces `_run` exactly (checked in gate GC2).
    """
    import core.amber as am
    from openmm import unit
    import openmm

    t0 = time.time()
    pos = np.asarray(pos_nm, float)
    ref = np.asarray(ref_nm, float)
    ctx = H.context
    ctx.setPositions(pos * unit.nanometer)
    for j, idx in enumerate(H._restraint_idx):
        H._rest_force.setParticleParameters(j, idx, ref[idx].tolist())
    H._rest_force.updateParametersInContext(ctx)
    e0 = am._energy_kcal(H)
    if steps >= 0:
        ctx.setParameter("k_rest", float(k_restraint) * am._K_SCALE)
        openmm.LocalEnergyMinimizer.minimize(ctx, float(tolerance), int(steps))
    ctx.setParameter("k_rest", 0.0)
    st = ctx.getState(getEnergy=True, getPositions=True)
    energy = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    out_nm = np.asarray(st.getPositions(asNumpy=True).value_in_unit(unit.nanometer), float)
    out = am._split(H, out_nm)
    out.update({"energy": float(energy), "energy_initial": float(e0),
                "k_restraint": float(k_restraint), "steps": int(steps),
                "tolerance": float(tolerance), "positions_nm": out_nm})
    d = (out_nm[H._restraint_idx] - ref[H._restraint_idx]) * 10.0
    out["restraint_rmsd"] = float(np.sqrt((d ** 2).sum(1).mean()))
    out["restraint_max"] = float(np.sqrt((d ** 2).sum(1)).max())
    out.update(am.convergence_flags(energy, e0))
    #: a wall-clock proxy for "the bound bit": a call that returns in far more than the typical
    #: time is the one that ran to `steps`.  Recorded, never used to change a number.
    out["wall"] = time.time() - t0
    out["hit_cap"] = bool(int(steps) > 0 and out["wall"] > 25.0)
    return out


def _record(r, coords_in, nat, seq):
    bb = {a: np.asarray(r["backbone"][a], float) for a in ATOMS if a in r["backbone"]}
    ca = np.asarray(r["ca"], float)
    return {"rmsd": float(I.ca_rmsd(ca, nat)),
            "input_rmsd": float(I.ca_rmsd(coords_in["CA"], nat)),
            "ca_disp": float(I.ca_rmsd(ca, coords_in["CA"])),
            "energy": float(r["energy"]), "energy_initial": float(r["energy_initial"]),
            "converged": bool(r["converged"]),
            "restraint_rmsd": float(r["restraint_rmsd"]),
            "restraint_max": float(r["restraint_max"]),
            "valid": EL.panel(bb, seq),
            "hit_cap": bool(r.get("hit_cap", False)),
            "wall": round(float(r["wall"]), 2)}


#: The rotated-frame null is EXACTLY ZERO by construction, so it is a numerical-floor
#: measurement, not an effect estimate: it is run on the first `N_FRAME_NULL` targets and
#: reported with its MAXIMUM (BRIEF section 7), never its mean.  Recorded as a compute choice.
N_FRAME_NULL = 15


#: The two softest rungs of the PRE-REGISTERED ladder, run UNCAPPED (`steps = 0`, the deployed
#: protocol) in a second pass.  They were dropped from pass 1 on a timing claim that I later
#: REFUTED myself (see agentC_FINDINGS.md section 4.8): the apparent non-termination was CPU
#: starvation on a shared box, not a property of the minimiser.  Restoring them restores the
#: pre-registration.
SOFT_LADDER = (1.0, 3.0)


def run_soft_target(t):
    import torsion_lib2 as tl2
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    nat = np.asarray(u["nat_ca"], float)
    avg, _C, _dev = averaged_backbone_from(np.asarray(W, float), np.asarray(PHI, float),
                                           np.asarray(PSI, float))
    coords = {a: np.asarray(avg[a], float) for a in ATOMS if a in avg}
    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    H = builder(seq, rep, SETS["bb"])
    pos0 = H._assemble(H._heavy_positions(coords, chi1=None))
    arms = {}
    for k in SOFT_LADDER:
        r = _run_ref(H, pos0, pos0, k, steps=0)      # UNCAPPED, the deployed protocol
        arms[f"k{int(k)}"] = _record(r, coords, nat, seq)
    drop_builders()
    return {"pdb": pdb, "n": n, "fold": int(t["fold"]), "seq": seq,
            "input_rmsd": float(I.ca_rmsd(coords["CA"], nat)),
            "input_valid": EL.panel(coords, seq), "arms": arms}


def run_soft(out="agentC_pareto_soft.json", verbose=True):
    import json
    tg = I.targets()
    path = os.path.join(CL.RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            rows = json.load(open(path)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    cfg = {"SOFT_LADDER": list(SOFT_LADDER), "STEPS": 0, "TOL": TOL,
           "note": "pass 2: the two softest PRE-REGISTERED rungs, uncapped"}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        PL.mem_ok(0.55)
        rows.append(run_soft_target(t))
        if verbose and len(rows) % 10 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
        _write(out, rows, cfg, len(tg))
    _write(out, rows, cfg, len(tg))
    return rows


def run_target(t, do_frame_null=True):
    import torsion_lib2 as tl2

    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    nat = np.asarray(u["nat_ca"], float)
    avg, _C, _dev = averaged_backbone_from(np.asarray(W, float), np.asarray(PHI, float),
                                           np.asarray(PSI, float))
    coords = {a: np.asarray(avg[a], float) for a in ATOMS if a in avg}

    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    Hb = builder(seq, rep, SETS["bb"])
    heavy = Hb._heavy_positions(coords, chi1=None)
    pos0 = Hb._assemble(heavy)

    arms = {}
    pos_of = {}
    for k in K_LADDER:
        r = _run_ref(Hb, pos0, pos0, k)
        arms[f"k{int(k)}"] = _record(r, coords, nat, seq)
        pos_of[f"k{int(k)}"] = r["positions_nm"]

    #: STAGED -- continue the k = 300 structure at k = 30, restrained to ITSELF (a genuine second
    #: stage: the fold is first frozen while strain relaxes, then softened).
    r = _run_ref(Hb, pos_of["k300"], pos_of["k300"], 30.0)
    arms["stage_300_30"] = _record(r, coords, nat, seq)

    #: PULL-BACK -- take the k = 30 structure and re-minimise it hard against the ORIGINAL input.
    #: A validity-preserving projection that stays at an ff14SB minimum instead of interpolating.
    r = _run_ref(Hb, pos_of["k30"], pos0, 1000.0)
    arms["pullback_30_1000"] = _record(r, coords, nat, seq)

    #: RESTRAINED-SET VARIANTS.  Only the restraint's atom list changes.
    for nm, aset in (("caonly", SETS["caonly"]), ("bbo", SETS["bbo"])):
        H2 = builder(seq, rep, aset)
        p2 = H2._assemble(H2._heavy_positions(coords, chi1=None))
        ks = (30.0, 100.0, 300.0) if nm == "caonly" else (30.0,)   # the full declared set
        for k in ks:
            r = _run_ref(H2, p2, p2, k)
            arms[f"{nm}_k{int(k)}"] = _record(r, coords, nat, seq)

    #: BLEND -- a free, AMBER-free step-back of the k = 30 output toward its input.  Declared in
    #: the pre-registration as a projection CANDIDATE expected to fail the validity axis.
    from core import amber as am
    k30_bb = {a: np.asarray(arms["k30"], float) for a in ()}    # placeholder, replaced below
    r30_bb = am._split(Hb, pos_of["k30"])["backbone"]
    for al in BLEND:
        bl = {a: (1.0 - al) * np.asarray(r30_bb[a], float) + al * coords[a]
              for a in r30_bb if a in coords}
        for a in coords:
            bl.setdefault(a, coords[a])
        arms[f"blend{int(al*100)}"] = {
            "rmsd": float(I.ca_rmsd(bl["CA"], nat)),
            "input_rmsd": float(I.ca_rmsd(coords["CA"], nat)),
            "ca_disp": float(I.ca_rmsd(bl["CA"], coords["CA"])),
            "energy": arms["k30"]["energy"], "energy_initial": arms["k30"]["energy_initial"],
            "converged": arms["k30"]["converged"],
            "restraint_rmsd": float("nan"), "restraint_max": float("nan"),
            "valid": EL.panel({a: bl[a] for a in ATOMS if a in bl}, seq), "wall": 0.0}

    rec = {"pdb": pdb, "n": n, "fold": int(t["fold"]), "seq": seq,
           "input_rmsd": float(I.ca_rmsd(coords["CA"], nat)),
           "input_valid": EL.panel(coords, seq), "arms": arms}

    #: THE ROTATED-FRAME NULL -- exactly zero by construction (ff14SB, the restraint, the
    #: projection and `ca_rmsd` are all rigid-invariant).  Reported with its MAXIMUM.
    if do_frame_null:
        rng = SD.stable_rng(pdb, "s19C_frame")
        Q, tvec = random_rigid(rng)
        cr = {a: (np.asarray(coords[a], float) @ Q.T) + tvec for a in coords}
        p3 = Hb._assemble(Hb._heavy_positions(cr, chi1=None))
        r = _run_ref(Hb, p3, p3, K_INCUMBENT)
        rec["frame_null"] = {"rmsd": float(I.ca_rmsd(np.asarray(r["ca"], float), nat)),
                             "energy": float(r["energy"]),
                             "d_rmsd": float(I.ca_rmsd(np.asarray(r["ca"], float), nat)
                                             - arms["k30"]["rmsd"]),
                             "d_energy": float(r["energy"] - arms["k30"]["energy"])}
    drop_builders()
    return rec


def run(targets=None, out="agentC_pareto.json", verbose=True):
    import json
    tg = targets if targets is not None else I.targets()
    path = os.path.join(CL.RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            rows = json.load(open(path)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    cfg = {"K_LADDER": list(K_LADDER), "STEPS": STEPS, "TOL": TOL, "BLEND": list(BLEND),
           "SETS": {k: list(v) for k, v in SETS.items()}, "K_INCUMBENT": K_INCUMBENT}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        PL.mem_ok(0.55)
        rows.append(run_target(t, do_frame_null=(len(rows) < N_FRAME_NULL)))
        if verbose:
            a = rows[-1]["arms"]
            print(f"  {len(rows)}/{len(tg)} {t['pdb']} in={rows[-1]['input_rmsd']:.3f} "
                  f"k30={a['k30']['rmsd']:.3f} k300={a['k300']['rmsd']:.3f} "
                  f"ca30={a['caonly_k30']['rmsd']:.3f} pb={a['pullback_30_1000']['rmsd']:.3f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        _write(out, rows, cfg, len(tg))
    _write(out, rows, cfg, len(tg))
    return rows


def _write(out, rows, cfg, n_expected):
    import json
    fn = [r["frame_null"]["d_rmsd"] for r in rows if "frame_null" in r]
    fe = [r["frame_null"]["d_energy"] for r in rows if "frame_null" in r]
    obj = {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg),
           "n_rows": len(rows), "n_expected": int(n_expected),
           "complete": bool(len(rows) >= int(n_expected))}
    if fn:
        obj["frame_null"] = {"n": len(fn),
                             "max_abs": float(np.max(np.abs(fn))),
                             "mean_abs": float(np.mean(np.abs(fn))),
                             "max_abs_energy": float(np.max(np.abs(fe)))}
    p = os.path.join(CL.RESULTS, out if out.endswith(".json") else out + ".json")
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, p)
    return p


def gate_GC3(n_targets=5, cap=None, verbose=True):
    """GC3 -- the declared iteration bound is INERT for the incumbent protocol.

    Capped k = 30 must reproduce uncapped k = 30 bit-for-bit on every checked target, or the
    bound is not a compute decision but a change of operator.
    """
    import torsion_lib2 as tl2
    from core import amber as am
    PL.mem_ok(0.55)
    cap = STEPS if cap is None else int(cap)
    worst_ca = worst_e = 0.0
    rows = []
    for t in I.targets()[:n_targets]:
        seq = t["seq"]
        W, PHI, PSI, u = top75_windows(t["pdb"])
        avg, _C, _d = averaged_backbone_from(np.asarray(W, float), np.asarray(PHI, float),
                                             np.asarray(PSI, float))
        coords = {a: np.asarray(avg[a], float) for a in ATOMS if a in avg}
        tab = tl2.library_for(seq, 4, seq)
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        H = builder(seq, rep, SETS["bb"])
        pos = H._assemble(H._heavy_positions(coords, chi1=None))
        a = _run_ref(H, pos, pos, K_INCUMBENT, steps=0, tolerance=TOL)
        b = _run_ref(H, pos, pos, K_INCUMBENT, steps=cap, tolerance=TOL)
        dca = float(np.abs(np.asarray(a["ca"], float) - np.asarray(b["ca"], float)).max())
        de = abs(float(a["energy"]) - float(b["energy"]))
        worst_ca = max(worst_ca, dca); worst_e = max(worst_e, de)
        rows.append({"pdb": t["pdb"], "d_ca_A": dca, "d_energy_kcal": de,
                     "wall_uncapped": round(a["wall"], 2), "wall_capped": round(b["wall"], 2)})
        drop_builders()
    out = {"cap": int(cap), "n_targets": int(n_targets), "rows": rows,
           "max_d_ca_A": worst_ca, "max_d_energy_kcal": worst_e,
           "passed": bool(worst_ca < 1e-9 and worst_e < 1e-9)}
    if verbose:
        print(f"GC3  cap={cap}: max |delta Ca| {worst_ca:.3e} A, max |delta E| {worst_e:.3e} "
              f"kcal/mol over {n_targets} targets -> {'PASS' if out['passed'] else 'FAIL'}")
        for r in rows:
            print(f"     {r['pdb']}  uncapped {r['wall_uncapped']:>6.2f}s  "
                  f"capped {r['wall_capped']:>6.2f}s  dCa {r['d_ca_A']:.2e}")
    return out


def gate_GC2(verbose=True):
    """GC2 -- `_run_ref(H, pos, pos, k)` IS `core.amber._run(H, heavy, k)`, and the patched
    builder really changes the restrained set without changing anything else."""
    import torsion_lib2 as tl2
    from core import amber as am
    PL.mem_ok(0.55)
    t = I.targets()[0]
    seq = t["seq"]
    W, PHI, PSI, u = top75_windows(t["pdb"])
    avg, _C, _d = averaged_backbone_from(np.asarray(W, float), np.asarray(PHI, float),
                                         np.asarray(PSI, float))
    coords = {a: np.asarray(avg[a], float) for a in ATOMS if a in avg}
    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    H = builder(seq, rep, SETS["bb"])
    heavy = H._heavy_positions(coords, chi1=None)
    a = am._run(H, heavy, K_INCUMBENT, STEPS, TOL, False)
    b = _run_ref(H, H._assemble(heavy), H._assemble(heavy), K_INCUMBENT)
    dca = float(np.abs(np.asarray(a["ca"], float) - np.asarray(b["ca"], float)).max())
    de = abs(float(a["energy"]) - float(b["energy"]))
    #: and the deployed path, which is what every Sprint-15..18 AMBER number was made with.
    #: `refine_coords` goes through `core.amber.builder_for`, whose `memory_guard()` refuses to
    #: build a Context above 92% physical memory; on a shared box that leg can be unavailable.
    #: It is then recorded as SKIPPED rather than silently passed.
    try:
        c = am.refine_coords(seq, rep, coords, k_restraint=K_INCUMBENT, steps=STEPS,
                             tolerance=TOL, threads=1, memo=False)
        dca2 = float(np.abs(np.asarray(a["ca"], float) - np.asarray(c["ca"], float)).max())
    except MemoryError as e:
        dca2 = float("nan")
        print(f"GC2  deployed-path leg SKIPPED: {e}")
    H2 = builder(seq, rep, SETS["caonly"])
    nres = len(seq)
    ok_sets = (len(H._restraint_idx) == 3 * nres and len(H2._restraint_idx) == nres
               and am.RESTRAINED_BACKBONE == ("N", "CA", "C"))
    out = {"max_dca_runref_vs_run_A": dca, "d_energy_kcal": de,
           "max_dca_run_vs_refine_coords_A": dca2,
           "restraint_sizes": [len(H._restraint_idx), len(H2._restraint_idx), 3 * nres],
           "module_constant_restored": bool(am.RESTRAINED_BACKBONE == ("N", "CA", "C")),
           "deployed_leg_run": bool(np.isfinite(dca2)),
           "passed": bool(dca < 1e-9 and de < 1e-9 and ok_sets
                          and (not np.isfinite(dca2) or dca2 < 1e-9))}
    drop_builders()
    if verbose:
        print(f"GC2  _run_ref vs _run {dca:.3e} A / {de:.3e} kcal | vs refine_coords {dca2:.3e} A"
              f" | restraint sizes {out['restraint_sizes']} | "
              f"{'PASS' if out['passed'] else 'FAIL'}")
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--gate", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--soft", action="store_true")
    ap.add_argument("--out", default="agentC_pareto.json")
    a = ap.parse_args()
    if a.soft:
        run_soft()
        sys.exit(0)
    if a.gate:
        g3 = gate_GC3()
        g2 = gate_GC2()
        import json as _j
        with open(os.path.join(CL.RESULTS, "agentC_gates_amber.json"), "w") as fh:
            _j.dump({"GC2": g2, "GC3": g3, "complete": True,
                     "passed": bool(g2["passed"] and g3["passed"])}, fh, indent=1)
        sys.exit(0)
    tg = I.targets()
    out = a.out
    if a.smoke:
        tg, out = tg[:4], "_SMOKE_agentC_pareto.json"
    elif a.limit:
        tg, out = tg[:a.limit], "_PARTIAL_" + a.out
    run(tg, out=out)
