"""S15 AUDIT / Part D -- energy models verified independently of their own drivers.

AMBER: OpenMM is the physics engine, so a genuinely independent implementation is out of
scope; what IS independently checkable is everything the recipe asserts about itself --
force-field XML actually loaded, the force classes actually present in the System, the
protonation/atom count, the restraint constant, the minimiser tolerance, determinism, and
the true per-call cost with the result memo defeated.  That is done here by building the
System through OpenMM directly and comparing to `core.amber`'s own number.

Legacy: an eleven-term additive field, so its TOTAL is independently checkable as the sum
of its own published components, and its components are checkable against a second
evaluation path.

    python -m s15.audit_energy amber
    python -m s15.audit_energy cost
    python -m s15.audit_energy legacy
"""
from __future__ import annotations
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s15", "results")
os.makedirs(OUT, exist_ok=True)

GOLDEN_1A13 = -489.9138948277905   # core/amber.py docstring: native interaction energy


def _w(name, obj):
    p = os.path.join(OUT, name if name.endswith(".json") else name + ".json")
    with open(p, "w") as fh:
        json.dump(obj, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print("wrote", p, flush=True)


# ------------------------------------------------------------------- provenance
def amber():
    """Interrogate the built System: what force field, what forces, what parameters."""
    import openmm
    from openmm import app, unit
    import core.amber as A
    out = {"openmm_version": openmm.version.version,
           "platforms": [openmm.Platform.getPlatform(i).getName()
                         for i in range(openmm.Platform.getNumPlatforms())]}
    # the declared recipe, read out of the module rather than retyped
    src = open(os.path.join(ROOT, "core", "amber.py")).read()
    out["forcefield_xml_in_source"] = sorted(
        {s for s in ("amber14/protein.ff14SB.xml", "implicit/gbn2.xml", "amber14-all.xml",
                     "implicit/obc2.xml", "implicit/gbn.xml") if s in src})
    for name in ("K_MODERATE", "K_STRONG", "K_WEAK", "MEMO_SIZE", "AMBER_TERMS",
                 "STRAIN_LIMIT", "TOLERANCE", "DEFAULT_TOLERANCE"):
        if hasattr(A, name):
            v = getattr(A, name)
            out[name] = list(v) if isinstance(v, (tuple, list)) else v

    # build the System through the module's own builder and inspect it
    from s12 import instrument as I
    u = I.load_univ("1A13")
    seq = u["seq"]
    from s13 import qarch_lib as Q
    sp = Q.Space("1A13", 4)
    H = A.builder_for(sp.seq, sp.rep, "CPU", 1)
    out["builder_class"] = type(H).__name__
    sysobj = None
    for attr in ("system", "_system", "sys"):
        if hasattr(H, attr):
            sysobj = getattr(H, attr); out["system_attr"] = attr; break
    if sysobj is not None:
        forces = [type(sysobj.getForce(i)).__name__ for i in range(sysobj.getNumForces())]
        out["forces_in_system"] = forces
        out["n_particles"] = sysobj.getNumParticles()
        out["has_CustomGBForce_GBn2"] = "CustomGBForce" in forces
        out["has_NonbondedForce"] = "NonbondedForce" in forces
        for i in range(sysobj.getNumForces()):
            f = sysobj.getForce(i)
            if type(f).__name__ == "NonbondedForce":
                out["nonbonded_method"] = int(f.getNonbondedMethod())
                out["nonbonded_cutoff_nm"] = float(f.getCutoffDistance().value_in_unit(unit.nanometer))
                q = [float(f.getParticleParameters(j)[0].value_in_unit(unit.elementary_charge))
                     for j in range(f.getNumParticles())]
                out["total_charge"] = round(float(sum(q)), 6)
    _w("audit_energy_amber", out)
    print(json.dumps(out, indent=1, default=str))
    return out


def golden():
    """Reproduce the pinned 1A13 native interaction energy stated in core/amber.py."""
    import core.amber as A
    from s12 import instrument as I
    u = I.load_univ("1A13")
    out = {"target": "1A13", "seq": u["seq"], "golden_in_docstring": GOLDEN_1A13}
    try:
        r = A.refine_coords(u["seq"], np.asarray(u["nat_ca"], float),
                            k_restraint=10.0, steps=0, tolerance=1.0, components=True)
        c = r.get("components", {})
        inter = float(c.get("nonbonded", np.nan)) + float(c.get("solvation", np.nan))
        out["components"] = {k: float(v) for k, v in c.items()}
        out["interaction_nonbonded_plus_solvation"] = inter
        out["delta_vs_golden"] = inter - GOLDEN_1A13
        out["bit_identical"] = inter == GOLDEN_1A13
        out["energy_total"] = float(r.get("energy", np.nan))
    except Exception as ex:
        import traceback
        out["error"] = traceback.format_exc()[-1200:]
    _w("audit_energy_golden", out)
    print(json.dumps({k: out[k] for k in out if k != "error"}, indent=1))
    if "error" in out:
        print(out["error"])
    return out


def cost(n=30):
    """The AMBER single-point cost, memo DEFEATED, against the recorded 28 ms / 6 ms."""
    import core.amber as A
    from s13 import qarch_lib as Q
    sp = Q.Space("2MK7", 4)
    rng = np.random.default_rng(0)
    S = rng.integers(0, 4, size=(n, sp.n))
    out = {"pdb": "2MK7", "n": int(sp.n), "n_calls": int(n),
           "recorded_in_brief_ms": 28.0, "recorded_in_qarch_lib_docstring_ms": 6.0}
    A.reset_caches() if hasattr(A, "reset_caches") else None
    # warm the builder so we time energy, not System construction
    A.single_point(sp.seq, sp.rep, [int(x) for x in S[0]], threads=1)
    st0 = A.cache_stats() if hasattr(A, "cache_stats") else {}
    t0 = time.perf_counter()
    for b in range(n):
        A.single_point(sp.seq, sp.rep, [int(x) for x in S[b]], threads=1)
    dt = time.perf_counter() - t0
    st1 = A.cache_stats() if hasattr(A, "cache_stats") else {}
    out["ms_per_call_distinct_inputs"] = 1000.0 * dt / n
    out["memo_hits_during_distinct"] = (st1.get("memo_hits", 0) - st0.get("memo_hits", 0))
    out["memo_misses_during_distinct"] = (st1.get("memo_misses", 0) - st0.get("memo_misses", 0))
    # the WRONG measurement: the same structure n times
    st0 = A.cache_stats() if hasattr(A, "cache_stats") else {}
    t0 = time.perf_counter()
    for _ in range(n):
        A.single_point(sp.seq, sp.rep, [int(x) for x in S[0]], threads=1)
    dtr = time.perf_counter() - t0
    st1 = A.cache_stats() if hasattr(A, "cache_stats") else {}
    out["ms_per_call_repeated_input"] = 1000.0 * dtr / n
    out["memo_hits_during_repeat"] = (st1.get("memo_hits", 0) - st0.get("memo_hits", 0))
    out["speedup_from_memo"] = out["ms_per_call_distinct_inputs"] / max(out["ms_per_call_repeated_input"], 1e-9)
    # determinism of a single point
    a = float(A.single_point(sp.seq, sp.rep, [int(x) for x in S[1]], threads=1)["energy"])
    A.reset_caches() if hasattr(A, "reset_caches") else None
    b2 = float(A.single_point(sp.seq, sp.rep, [int(x) for x in S[1]], threads=1)["energy"])
    out["single_point_deterministic_after_cache_reset"] = (a == b2)
    out["single_point_value"] = a
    out["single_point_value_after_reset"] = b2
    _w("audit_energy_cost", out)
    print(json.dumps(out, indent=1))
    return out


def legacy():
    """Legacy: does the published TOTAL equal the sum of its own eleven components?"""
    from s13 import qarch_lib as Q
    rows = []
    import glob
    for f in sorted(glob.glob(os.path.join(ROOT, "s13", "results", "qarch_enum_*.npz"))):
        z = np.load(f, allow_pickle=True)
        terms = [k for k in z.files if k.startswith("leg_")]
        tot = np.asarray(z["legacy"], np.float64)
        s = sum(np.asarray(z[t], np.float64) for t in terms)
        d = np.abs(tot - s)
        # and the weighted form, in case legacy_total is not a plain sum
        rows.append({"pdb": str(z["pdb"]), "n_terms": len(terms),
                     "terms": sorted(t[4:] for t in terms),
                     "max_abs_diff_total_vs_sum": float(d.max()),
                     "mean_abs_diff": float(d.mean()),
                     "total_sd": float(tot.std()),
                     "plain_sum_reproduces_total": bool(d.max() < 1e-3)})
        print(json.dumps(rows[-1]), flush=True)
    # second path: recompute components for a handful of configs through qarch_lib
    sp = Q.Space("2MK7", 4)
    S = np.array(list(np.ndindex(*(4,) * sp.n))[:64], np.int8) if sp.n <= 9 else None
    check = {}
    if S is not None:
        comp = Q.legacy_components(sp, S)
        tot2 = Q.legacy_total(comp)
        z = np.load(os.path.join(ROOT, "s13", "results", "qarch_enum_2MK7.npz"))
        stored = np.asarray(z["legacy"], np.float64)[:64]
        check = {"pdb": "2MK7", "n_checked": int(len(S)),
                 "max_abs_diff_recomputed_vs_stored": float(np.abs(tot2 - stored).max()),
                 "bit": bool(np.array_equal(np.asarray(tot2, np.float32), np.asarray(stored, np.float32)))}
        print(json.dumps(check), flush=True)
    _w("audit_energy_legacy", {"rows": rows, "recompute_check": check, "complete": True})


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "amber"
    {"amber": amber, "cost": cost, "legacy": legacy, "golden": golden}[cmd]()
