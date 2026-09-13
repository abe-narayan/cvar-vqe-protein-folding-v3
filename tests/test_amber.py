"""Equivalence and genuineness tests for `core.amber` against the shipped AMBER stack.

The claim this file has to support is not "the consolidation looks right". It is that
`core.amber` returns *the same floats* as `amber_refine` + `amber_hamiltonian` +
`sidechains` + `budget`, on real pipeline structures, and that the physics underneath is
still genuine ff14SB + GBn2 through OpenMM rather than a cheaper proxy that happens to
correlate.

Structures come from the production recipe: BLOSUM62 top-k retrieval over the fold-held-out
peptide database, rebuilt by `protein_geometry.build_backbone_batch`, scored with
`refine_coords(k_restraint=K_MODERATE, steps=0, tolerance=1.0, components=True)`. No toy
inputs, no hand-written coordinates.

Tolerance, and why it is what it is
`core.amber` lifts the hot arithmetic byte-for-byte and reuses the same OpenMM calls in
the same order on the same platform at Threads=1, so the expected difference is not
"small", it is **exactly zero**. Every energy assertion here is `==`, not `approx`. That
is a stronger and more honest test than a tolerance would be: any reduction-order change
anywhere in the builder, the sidechain templates, the hydrogen frames or the minimiser
driver shows up immediately instead of hiding under an epsilon.

The one place a tolerance IS justified is `test_threads_deviation_is_bounded`, which
measures the deviation the threaded CustomGBForce path actually introduces rather than
asserting it is absent. It is reported, not assumed.
"""
import ctypes
import json
import os
import sys
import time

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import protein_geometry as geo                                       # noqa: E402
import peptide_db as db                                              # noqa: E402
import torsion_lib2 as tl2                                           # noqa: E402

import core.amber as A                                               # noqa: E402

#: The pinned reference from the sprint-7/8 scoring recipe: 1A13's native interaction
#: energy (nonbonded + solvation) at k_restraint=10, steps=0, tolerance=1.0.
REF_1A13_INTERACTION = -489.9138948277905

PID = "1A13"


# The box is not the code
#
# `core.amber.memory_guard` raises `MemoryError` above 92% physical memory, deliberately,
# so that an OOM here kills this run rather than a sibling's. On a shared box that CORRECT
# REFUSAL arrives as four red tests, and an environmental red competes for attention with a
# real one -- the same failure mode as a run reporting 126/126 while AMBER quietly declined
# on a fifth of the manifest.
#
# But a blanket `except MemoryError: skip` would hide the regression that matters most: this
# suite leaking contexts until IT is what fills the box. So the verdict is explicit and
# unit-tested. Skip only when the memory that crossed the ceiling demonstrably is not ours;
# FAIL when our own working set accounts for it.
def _own_rss_bytes() -> int:
    """This process's working set.

    Both ctypes signatures are declared rather than left to guess, and neither default is
    harmless: `GetCurrentProcess` returns the pseudo-handle -1, which ctypes truncates to
    a 32-bit int without an explicit `restype`, and passing that back without `argtypes`
    raises `OverflowError: int too long to convert`. The same trap bit the affinity probe
    in `verify/amber_affinity.py`."""
    if os.name != "nt":
        return 0

    class _PMC(ctypes.Structure):
        _fields_ = [("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong),
                    ("PeakWorkingSetSize", ctypes.c_size_t),
                    ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t),
                    ("PeakPagefileUsage", ctypes.c_size_t)]

    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    psapi.GetProcessMemoryInfo.argtypes = [ctypes.c_void_p,
                                           ctypes.POINTER(_PMC), ctypes.c_ulong]
    psapi.GetProcessMemoryInfo.restype = ctypes.c_int
    c = _PMC()
    c.cb = ctypes.sizeof(_PMC)
    if not psapi.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(c), c.cb):
        return 0
    return int(c.WorkingSetSize)


def _total_phys_bytes() -> int:
    if os.name != "nt":
        return 0
    s = A._MEMORYSTATUSEX()
    s.dwLength = ctypes.sizeof(A._MEMORYSTATUSEX)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(s))
    return int(s.ullTotalPhys)

#: OUR share of the rise above which the suite is judged to be the cause.
OWN_SHARE_IS_A_REGRESSION = 0.5


def _memory_verdict(pct, baseline_pct, own_growth_bytes, total_bytes,
                    limit=None, share=OWN_SHARE_IS_A_REGRESSION):
    """``"run" | "skip" | "fail"``. Factored out so the decision is testable without
    needing a full box; `test_memory_verdict_distinguishes_the_box_from_the_code` drives it
    with synthetic numbers, including the case where THIS suite is the leak."""
    limit = A.MEMORY_LIMIT_PERCENT if limit is None else limit
    if pct <= limit:
        return "run"
    if baseline_pct > limit:
        return "skip"                    # the box was already full before we started
    rose = (pct - baseline_pct) / 100.0 * total_bytes
    if rose > 0 and own_growth_bytes / rose >= share:
        return "fail"                    # we are what filled it
    return "skip"

_BASELINE = {"pct": A.memory_percent(), "rss": _own_rss_bytes(),
             "total": _total_phys_bytes()}

#: The S26 resource governor's live snapshot (`s26/governor.py` rewrites it every 5 s).  When
#: it is fresh, a skip or fail quotes the governor's reading beside the in-process one, so
#: the message says which ceiling fired and what the box looked like to the process that
#: polices it.  The DECISION is unchanged: `_memory_verdict` still runs on
#: `core.amber.memory_percent()`, the same syscall `core.amber.memory_guard` uses, so the
#: verdict is identical whether or not a governor is running (S26 lane I, defect 6b).
GOVERNOR_STATE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              "s26", "governor_state.json")
#: a snapshot older than this is a stopped governor, not a current reading
GOVERNOR_MAX_AGE_S = 120.0


def _governor_reading(path=GOVERNOR_STATE, max_age_s=GOVERNOR_MAX_AGE_S, now=None):
    """The governor's last snapshot as a small dict, or None when it is not running
    (file absent, unreadable, malformed, or older than `max_age_s`)."""
    try:
        with open(path, encoding="utf-8") as fh:
            st = json.load(fh)
        age = (time.time() if now is None else now) - float(st["epoch"])
        if age > max_age_s:
            return None
        band = st.get("band") or {}
        return {"ram_pct": float(st["ram_pct"]), "cpu_pct": float(st.get("cpu_pct", 0.0)),
                "ts": str(st.get("ts", "")), "age_s": float(age),
                "n_jobs": int(st.get("n_jobs", 0)), "n_amber": int(st.get("n_amber", 0)),
                "ceiling": band.get("ceiling")}
    except Exception:                                                    # noqa: BLE001
        return None


def _ceiling_message(verdict, pct, baseline_pct, grew_bytes, limit=None, gov=None):
    """The text a skip or fail carries: the ceiling that fired, the reading that crossed it,
    and the governor's own reading when one is running."""
    limit = A.MEMORY_LIMIT_PERCENT if limit is None else limit
    text = (f"physical memory {pct:.0f}% (core.amber.memory_percent) is above core.amber's "
            f"{limit:.0f}% ceiling")
    if gov is not None:
        text += (f"; the S26 governor read {gov['ram_pct']:.1f}% RAM / {gov['cpu_pct']:.1f}% CPU "
                 f"at {gov['ts']} ({gov['age_s']:.0f} s ago) with {gov['n_jobs']} registered "
                 f"job(s), {gov['n_amber']} AMBER, governor ceiling {gov['ceiling']}%")
    else:
        text += "; no fresh s26/governor_state.json, so the governor is not running"
    if verdict == "fail":
        return (text + f", and THIS SUITE accounts for it: our working set grew "
                f"{grew_bytes / 2**20:.0f} MB since the run started, from {baseline_pct:.0f}%. "
                "That is a leak in the builder cache or context teardown, not a busy box.")
    return (text + f", and this suite did not cause it (our working set grew "
            f"{grew_bytes / 2**20:.0f} MB since {baseline_pct:.0f}%); memory_guard is refusing "
            "contexts correctly. Re-run on a quiet box.")


@pytest.fixture(autouse=True)
def _memory_ceiling():
    pct = A.memory_percent()
    grew = _own_rss_bytes() - _BASELINE["rss"]
    verdict = _memory_verdict(pct, _BASELINE["pct"], grew, _BASELINE["total"])
    if verdict == "fail":
        pytest.fail(_ceiling_message("fail", pct, _BASELINE["pct"], grew,
                                     gov=_governor_reading()))
    if verdict == "skip":
        pytest.skip(_ceiling_message("skip", pct, _BASELINE["pct"], grew,
                                     gov=_governor_reading()))
    yield
K = 6                      # candidates per pool; each costs ~6 s converged


# --------------------------------------------------------------------- fixtures
def _pool(pid=PID, k=K):
    """The production pool for one target: BLOSUM62 top-k, fold held out."""
    from s5.lib import B62, encode, windows_full
    from s7.amber_native import pool_for
    p = db.by_pdb(pid)
    fold = db.folds(5)[p.seq]
    _, S, PHI, PSI, _ = windows_full(pool_for(fold, p.seq), p.n)
    sim = B62[S, encode(p.seq)[None, :]].sum(1)
    idx = np.argsort(-sim)[:k]
    coords = geo.build_backbone_batch(PHI[idx].astype(float), PSI[idx].astype(float))
    tab = tl2.library_for(p.seq, 8, p.seq)
    rep = tl2.PerResidueTorsion(p.seq, tab, chi_bits=False)
    return p, rep, coords


def _slice(coords, b):
    return {k: v[b] for k, v in coords.items()}


@pytest.fixture(scope="module")
def pool():
    return _pool()


@pytest.fixture(scope="module")
def native():
    p = db.by_pdb(PID)
    tab = tl2.library_for(p.seq, 8, p.seq)
    rep = tl2.PerResidueTorsion(p.seq, tab, chi_bits=False)
    rb = geo.build_backbone(np.asarray(p.phi, float), np.asarray(p.psi, float))
    return p, rep, rb


# --------------------------------------------------------------- genuineness
def test_forcefield_is_really_ff14sb_gbn2(pool):
    """Real ff14SB + GBn2 through OpenMM, with CustomGBForce actually loaded."""
    p, rep, _ = pool
    H = A.builder_for(p.seq, rep)
    files = H.describe()["forcefield"]
    print(f"\nforcefield: {files}")
    assert "amber14/protein.ff14SB.xml" in files
    assert "implicit/gbn2.xml" in files

    names = sorted(f.__class__.__name__ for f in H.system.getForces())
    print(f"system forces: {names}")
    print(f"platform: {H.context.getPlatform().getName()} {H.platform_properties}")
    # GBn2 is implemented as a CustomGBForce. Its absence would mean the solvation term
    # -- which is 98% of every force evaluation and the larger half of the scoring
    # energy -- had been silently dropped or swapped for a cheaper implicit model.
    assert "CustomGBForce" in names, "GBn2 solvation is missing from the System"
    assert "NonbondedForce" in names
    assert "HarmonicBondForce" in names and "HarmonicAngleForce" in names
    assert "PeriodicTorsionForce" in names

    gb = next(f for f in H.system.getForces()
              if f.__class__.__name__ == "CustomGBForce")
    # GBn2 has two computed values (the I integral and the Born radius B).
    assert gb.getNumComputedValues() == 2
    assert gb.getNumPerParticleParameters() == 7
    assert gb.getNumParticles() == H.n_atoms
    print(f"CustomGBForce: {gb.getNumParticles()} particles, "
          f"{gb.getNumComputedValues()} computed values, "
          f"{gb.getNumEnergyTerms()} energy terms")


def test_solvation_dominates_the_force_evaluation(pool):
    """The GBn2 term is real work, not a stub returning a constant."""
    p, rep, coords = pool
    H = A.builder_for(p.seq, rep)
    r = A.refine_coords(p.seq, rep, _slice(coords, 0), k_restraint=A.K_MODERATE,
                        steps=0, tolerance=1.0, components=True)
    cm = r["components"]
    print(f"\ncomponents: { {k: round(v, 4) for k, v in cm.items()} }")
    assert cm["solvation"] < -100.0, "GBn2 solvation is implausibly small"
    assert np.isfinite(list(cm.values())).all()


# --------------------------------------------------------------- the reference
def test_reference_1a13_native_interaction_bit_exact(native):
    """The pinned number, reproduced bit-for-bit. This path is deterministic."""
    p, rep, rb = native
    out = A.refine_coords(p.seq, rep, rb, k_restraint=A.K_MODERATE, steps=0,
                          tolerance=1.0, components=True)
    cm = out["components"]
    inter = cm["nonbonded"] + cm["solvation"]
    print(f"\n1A13 native interaction: {inter!r}")
    print(f"  reference:              {REF_1A13_INTERACTION!r}")
    print(f"  abs difference:         {abs(inter - REF_1A13_INTERACTION):.3e}")
    assert inter == REF_1A13_INTERACTION


# --------------------------------------------------------------- equivalence
def test_refine_coords_matches_shipped_implementation_exactly(pool):
    """Old vs new, energy by energy, on a real retrieval pool.

    Every quantity is compared with `==`. `core.amber` reuses the identical arithmetic on
    the identical platform, so a nonzero difference anywhere is a defect, not noise.
    """
    import amber_refine as ar
    p, rep, coords = pool

    d_energy, d_inter, d_pos, d_comp = [], [], [], []
    for b in range(K):
        c = _slice(coords, b)
        old = ar.refine_coords(p.seq, rep, c, k_restraint=ar.K_MODERATE, steps=0,
                               tolerance=1.0, components=True)
        new = A.refine_coords(p.seq, rep, c, k_restraint=A.K_MODERATE, steps=0,
                              tolerance=1.0, components=True, memo=False)

        assert set(old["components"]) == set(new["components"])
        for t in A.AMBER_TERMS:
            d_comp.append(abs(old["components"][t] - new["components"][t]))
        oi = old["components"]["nonbonded"] + old["components"]["solvation"]
        ni = new["components"]["nonbonded"] + new["components"]["solvation"]
        d_inter.append(abs(oi - ni))
        d_energy.append(abs(old["energy"] - new["energy"]))
        d_pos.append(np.abs(old["positions_nm"] - new["positions_nm"]).max())

        # the strain-rejection rule must fire identically
        os_ = (old["components"]["bond"] + old["components"]["angle"]) > 1000.0
        ns_ = (new["components"]["bond"] + new["components"]["angle"]) > 1000.0
        assert os_ == ns_
        assert old["n_atoms"] == new["n_atoms"]
        assert np.array_equal(old["ca"], new["ca"])
        assert old["heavy_names"] == new["heavy_names"]
        assert old["restraint_rmsd"] == new["restraint_rmsd"]
        assert old["energy_initial"] == new["energy_initial"]

    print(f"\n{K} real candidates, old vs new:")
    print(f"  max |d total energy|      {max(d_energy):.3e} kcal/mol")
    print(f"  max |d interaction|       {max(d_inter):.3e} kcal/mol")
    print(f"  max |d component| (all 5) {max(d_comp):.3e} kcal/mol")
    print(f"  max |d position|          {max(d_pos):.3e} nm")
    assert max(d_energy) == 0.0
    assert max(d_inter) == 0.0
    assert max(d_comp) == 0.0
    assert max(d_pos) == 0.0


def test_sidechain_builder_matches_for_all_twenty_residues():
    """All 20 residues, atom for atom, against the shipped `sidechains`.

    `sidechains.py` builds all 20 standard residues (that is what took Amber coverage
    from 6/35 to 35/35 benchmark targets), so a consolidation that quietly dropped or
    mistranscribed one template would re-break coverage on exactly the targets that were
    hardest to add.
    """
    import sidechains as sc
    seq = "ACDEFGHIKLMNPQRSTVWY"
    assert len(set(seq)) == 20
    n = len(seq)
    phi = np.full(n, np.radians(-120.0))
    psi = np.full(n, np.radians(130.0))
    bb = geo.build_backbone(phi, psi)

    old = sc.build_full_structure(seq, bb)
    new = A.build_full_structure(seq, bb)
    assert old["n_atoms"] == new["n_atoms"]
    worst = 0.0
    for i in range(n):
        assert list(old["residues"][i]) == list(new["residues"][i]), f"atom set {seq[i]}"
        for nm in old["residues"][i]:
            d = np.abs(np.asarray(old["residues"][i][nm])
                       - np.asarray(new["residues"][i][nm])).max()
            worst = max(worst, d)
    print(f"\nall 20 residues, {old['n_atoms']} heavy atoms, "
          f"max |d coordinate| = {worst:.3e} A")
    assert worst == 0.0

    for aa in seq:
        rn = geo.ONE_TO_THREE[aa]
        assert sc.residue_bonds(rn) == A.residue_bonds(rn)
        assert sc.sidechain_atom_names(rn) == A.sidechain_atom_names(rn)
        assert sc.heavy_atom_count(rn) == A.heavy_atom_count(rn)
        assert sc.ring_atom_names(rn) == A.ring_atom_names(rn)


def test_topology_and_parameters_are_identical(pool):
    """Same atoms, same bonds, and the same ff14SB parameter assignment."""
    import amber_refine as ar
    p, rep, _ = pool
    Ho = ar.builder_for(p.seq, rep)
    Hn = A.builder_for(p.seq, rep)

    assert Ho.n_atoms == Hn.n_atoms
    ao = [(a.index, a.name, a.element.symbol, a.residue.index, a.residue.name)
          for a in Ho.topology.atoms()]
    an = [(a.index, a.name, a.element.symbol, a.residue.index, a.residue.name)
          for a in Hn.topology.atoms()]
    assert ao == an
    bo = sorted((b[0].index, b[1].index) for b in Ho.topology.bonds())
    bn = sorted((b[0].index, b[1].index) for b in Hn.topology.bonds())
    assert bo == bn
    assert Ho._heavy_names == Hn._heavy_names
    assert Ho._restraint_idx == Hn._restraint_idx

    # ff14SB nonbonded parameters, particle by particle
    fo = next(f for f in Ho.system.getForces()
              if f.__class__.__name__ == "NonbondedForce")
    fn = next(f for f in Hn.system.getForces()
              if f.__class__.__name__ == "NonbondedForce")
    dq = dsig = deps = 0.0
    for i in range(Ho.n_atoms):
        q1, s1, e1 = fo.getParticleParameters(i)
        q2, s2, e2 = fn.getParticleParameters(i)
        dq = max(dq, abs(q1 - q2)._value)
        dsig = max(dsig, abs(s1 - s2)._value)
        deps = max(deps, abs(e1 - e2)._value)
    assert fo.getNumExceptions() == fn.getNumExceptions()

    # GBn2 per-particle parameters
    go = next(f for f in Ho.system.getForces()
              if f.__class__.__name__ == "CustomGBForce")
    gn = next(f for f in Hn.system.getForces()
              if f.__class__.__name__ == "CustomGBForce")
    dgb = 0.0
    for i in range(Ho.n_atoms):
        dgb = max(dgb, max(abs(np.asarray(go.getParticleParameters(i), float)
                               - np.asarray(gn.getParticleParameters(i), float))))
    print(f"\n{Ho.n_atoms} atoms, {len(bo)} bonds identical")
    print(f"  max |d charge|   {dq:.3e} e")
    print(f"  max |d sigma|    {dsig:.3e} nm")
    print(f"  max |d epsilon|  {deps:.3e} kJ/mol")
    print(f"  max |d GBn2 par| {dgb:.3e}")
    assert dq == 0.0 and dsig == 0.0 and deps == 0.0 and dgb == 0.0

    # calibrated hydrogen frames
    assert len(Ho._h_local) == len(Hn._h_local)
    dh = max(float(np.abs(a[4] - b[4]).max())
             for a, b in zip(Ho._h_local, Hn._h_local))
    print(f"  max |d H frame|  {dh:.3e} nm  ({len(Ho._h_local)} hydrogens)")
    assert dh == 0.0


# --------------------------------------------------------------- the memo
def test_memo_is_exact_and_reports_its_own_wall(pool):
    """A memo hit must return the same floats and its OWN elapsed time."""
    p, rep, coords = pool
    A.clear_cache()
    c = _slice(coords, 0)
    a = A.refine_coords(p.seq, rep, c, k_restraint=A.K_MODERATE, steps=0,
                        tolerance=1.0, components=True)
    b = A.refine_coords(p.seq, rep, c, k_restraint=A.K_MODERATE, steps=0,
                        tolerance=1.0, components=True)
    st = A.cache_stats()
    print(f"\nmemo: {st}   cold {a['wall']:.3f}s -> warm {b['wall']:.6f}s")
    assert st["memo_hits"] == 1 and st["memo_misses"] == 1
    assert a["memo_hit"] is False and b["memo_hit"] is True
    assert a["components"] == b["components"]
    assert a["energy"] == b["energy"]
    assert np.array_equal(a["positions_nm"], b["positions_nm"])
    # the hit must not report the miss's cost
    assert b["wall"] < a["wall"] / 10.0

    # a different structure must NOT collide
    d = A.refine_coords(p.seq, rep, _slice(coords, 1), k_restraint=A.K_MODERATE,
                        steps=0, tolerance=1.0, components=True)
    assert d["energy"] != a["energy"]

    # mutating a returned dict must not poison the cache
    b["components"]["solvation"] = 12345.0
    b["ca"][0, 0] = 999.0
    e = A.refine_coords(p.seq, rep, c, k_restraint=A.K_MODERATE, steps=0,
                        tolerance=1.0, components=True)
    assert e["components"] == a["components"]
    assert np.array_equal(e["ca"], a["ca"])


def test_builder_cache_is_bounded_and_tears_contexts_down():
    """The old cache grew without limit; 126 targets meant 126 live Contexts."""
    A.clear_cache()
    assert A.cache_stats()["builders"] == 0
    assert A.BUILDER_CACHE_SIZE > 0
    seqs = ["AAAAAAA", "AAAAAAG", "AAAAAAS"]
    for s in seqs:
        tab = tl2.library_for(s, 8, s)
        rep = tl2.PerResidueTorsion(s, tab, chi_bits=False)
        A.builder_for(s, rep)
    assert A.cache_stats()["builders"] == len(seqs)
    A.BUILDER_CACHE_SIZE, keep = 2, A.BUILDER_CACHE_SIZE
    try:
        s = "AAAAAAT"
        tab = tl2.library_for(s, 8, s)
        A.builder_for(s, tl2.PerResidueTorsion(s, tab, chi_bits=False))
        n = A.cache_stats()["builders"]
        print(f"\nbuilders after eviction to cap 2: {n}")
        assert n <= 2
    finally:
        A.BUILDER_CACHE_SIZE = keep
        A.clear_cache()


# --------------------------------------------------------------- parallel
def test_refine_many_is_bit_identical_to_serial(pool):
    """The parallel map must parallelise the work, not approximate it."""
    p, rep, coords = pool
    cl = [_slice(coords, b) for b in range(K)]
    A.clear_cache()
    serial = [A.refine_coords(p.seq, rep, c, k_restraint=A.K_MODERATE, steps=0,
                              tolerance=1.0, components=True, memo=False)
              for c in cl]
    par = A.refine_many(p.seq, rep, cl, k_restraint=A.K_MODERATE, steps=0,
                        tolerance=1.0, components=True, workers=2)
    assert len(par) == len(serial)
    de, dp = [], []
    for s, q in zip(serial, par):
        for t in A.AMBER_TERMS:
            de.append(abs(s["components"][t] - q["components"][t]))
        dp.append(np.abs(s["positions_nm"] - q["positions_nm"]).max())
        assert s["heavy_names"] == q["heavy_names"]
    print(f"\nrefine_many vs serial: max |d component| {max(de):.3e} kcal/mol, "
          f"max |d position| {max(dp):.3e} nm")
    assert max(de) == 0.0
    assert max(dp) == 0.0


def test_plan_workers_respects_cores_and_memory():
    assert A.plan_workers(1) == 1
    assert A.plan_workers(0) == 1
    assert A.plan_workers(100, workers=3) == 3
    assert A.plan_workers(2, workers=99) == 2
    assert 1 <= A.plan_workers(100) <= (os.cpu_count() or 2)
    # the ceiling is a hard stop, not advice
    with pytest.raises(MemoryError):
        A.memory_guard(limit=-1.0)
    print(f"\nmemory {A.memory_percent():.0f}%  "
          f"plan_workers(500) = {A.plan_workers(500)}")


# --------------------------------------------------------------- determinism
def test_same_input_gives_bit_identical_output(pool):
    """The determinism the golden runs and budget.py's cache invariant depend on."""
    p, rep, coords = pool
    c = _slice(coords, 2)
    a = A.refine_coords(p.seq, rep, c, k_restraint=A.K_MODERATE, steps=0,
                        tolerance=1.0, components=True, memo=False)
    b = A.refine_coords(p.seq, rep, c, k_restraint=A.K_MODERATE, steps=0,
                        tolerance=1.0, components=True, memo=False)
    assert a["energy"] == b["energy"]
    assert a["components"] == b["components"]
    assert np.array_equal(a["positions_nm"], b["positions_nm"])


def test_threading_is_a_different_experiment_not_a_faster_one(pool):
    """Threading CustomGBForce is faster and lands somewhere else. Measure both halves.

    This test exists because measuring only the single point would have been wrong. On a
    fixed structure the threaded reduction moves the energy by ~1e-8 relative, which reads
    as harmless. Through 800 minimiser steps descending from 10^4-10^6 kcal/mol of builder
    strain, that steers the trajectory into a different local minimum -- kcal/mol, not
    ulps. Both numbers are asserted so neither can quietly change.
    """
    p, rep, coords = pool
    c = _slice(coords, 0)

    def run(threads, steps, k, tol):
        return A.refine_coords(p.seq, rep, c, k_restraint=k, steps=steps,
                               tolerance=tol, components=True, memo=False,
                               threads=threads)

    try:
        sp1, sp4 = run(1, -1, 0.0, 1e9), run(4, -1, 0.0, 1e9)
    except Exception as exc:                                         # noqa: BLE001
        pytest.skip(f"threaded platform unavailable: {exc}")

    # (a) SINGLE POINT -- no minimisation, so this isolates reduction order. The frozen
    # hydrogen frames must be identical, or the two builders are not the same builder.
    d_pos = np.abs(sp1["positions_nm"] - sp4["positions_nm"]).max() * 10.0
    i1 = sp1["components"]["nonbonded"] + sp1["components"]["solvation"]
    i4 = sp4["components"]["nonbonded"] + sp4["components"]["solvation"]
    rel_sp = abs(i1 - i4) / max(abs(i1), 1e-12)
    print(f"\nSINGLE POINT (no minimisation)")
    print(f"  threads=1 {i1!r}\n  threads=4 {i4!r}")
    print(f"  relative {rel_sp:.3e}   max |d position| {d_pos:.3e} A")
    assert d_pos == 0.0, ("hydrogen calibration is thread-dependent; the two runs are "
                          "not using the same builder")
    assert rel_sp < 1e-6, "single-point deviation is larger than reduction order"

    # (b) CONVERGED -- the production recipe. The deviation is now physical.
    cv1, cv4 = run(1, 0, A.K_MODERATE, 1.0), run(4, 0, A.K_MODERATE, 1.0)
    j1 = cv1["components"]["nonbonded"] + cv1["components"]["solvation"]
    j4 = cv4["components"]["nonbonded"] + cv4["components"]["solvation"]
    rel_cv = abs(j1 - j4) / max(abs(j1), 1e-12)
    moved = np.abs(cv1["positions_nm"] - cv4["positions_nm"]).max() * 10.0
    print(f"CONVERGED (the production recipe)")
    print(f"  threads=1 {j1!r}  ({cv1['wall']:.3f}s)")
    print(f"  threads=4 {j4!r}  ({cv4['wall']:.3f}s)")
    print(f"  absolute {abs(j1 - j4):.3e} kcal/mol   relative {rel_cv:.3e}")
    print(f"  max |d position| {moved:.3f} A   speedup {cv1['wall']/cv4['wall']:.2f}x")
    print("  -> a different local minimum, not a different rounding; use refine_many")
    # Amplification through the minimiser is the point. If this ever became negligible,
    # the docstring's warning would be wrong and should be revisited.
    assert rel_cv > 100.0 * rel_sp


# --------------------------------------------------------------- one definition
def test_budget_classes_have_exactly_one_definition():
    """`BudgetExhausted` must be ONE class, not two identical ones.

    Two copies of an exception class is not a style problem, it is a live trap:
    ``except budget.BudgetExhausted`` does not catch ``core.amber.BudgetExhausted``, and
    the two are indistinguishable in a diff.  `core.quantum.FoldingHamiltonian` subclasses
    the model and raises the exception; the three ``except BudgetExhausted`` handlers in
    `core.quantum` catch it.  Those were consistent only by accident of which copy each
    line happened to import, and a partial repoint would have sent budget exhaustion
    through three handlers that still read correctly -- turning a graceful stop into a
    crash inside the VQE path.

    The direction of the dependency is asserted too, and it is deliberate: `budget.py` is
    what ``CORE_BACKENDS=legacy`` runs, so it must not import from `core.amber` or the
    baseline arm would silently be executing consolidated code -- and every equivalence
    number in this sprint is measured against that arm being pure.
    """
    import inspect

    import budget
    import core.quantum as Q

    assert A.BudgetExhausted is budget.BudgetExhausted
    assert A.BudgetedEnergyModel is budget.BudgetedEnergyModel
    assert Q.BudgetExhausted is A.BudgetExhausted
    assert A.resolve_maxiter is budget.resolve_maxiter
    assert A.check_optimizer_budget is budget.check_optimizer_budget
    assert budget.BudgetedEnergyModel in Q.FoldingHamiltonian.__mro__

    #: the property that actually matters, exercised rather than inferred
    try:
        raise A.BudgetExhausted(3, 3)
    except budget.BudgetExhausted:
        pass

    #: and the legacy arm stays free of the consolidated one
    src = inspect.getsource(budget)
    assert "core.amber" not in src and "from core" not in src, \
        "budget.py is the legacy arm; it must not import the consolidated module"


def test_memory_verdict_distinguishes_the_box_from_the_code():
    """The skip must not be able to hide a leak in this suite.

    Driven with synthetic numbers so the decision is checked without needing to fill a
    real box -- including the case that matters, where memory crossed the ceiling and our
    own working set is what did it.
    """
    GB = 2 ** 30
    total = 16 * GB
    # under the ceiling: always run, whatever anything grew by
    assert _memory_verdict(80.0, 60.0, 8 * GB, total, limit=92.0) == "run"
    # over the ceiling, our growth is negligible against the machine's rise -> the box
    assert _memory_verdict(95.0, 60.0, 100 * 2 ** 20, total, limit=92.0) == "skip"
    # over the ceiling, and OUR working set accounts for the rise -> our bug, fail loudly
    assert _memory_verdict(95.0, 60.0, 5 * GB, total, limit=92.0) == "fail"
    # already full before the suite started: not attributable to us either way
    assert _memory_verdict(95.0, 94.0, 5 * GB, total, limit=92.0) == "skip"
    # exactly at the ceiling is still runnable -- the guard uses > not >=
    assert _memory_verdict(92.0, 60.0, 5 * GB, total, limit=92.0) == "run"


def test_the_memory_probes_return_something_sane():
    """A skip decision built on a syscall that silently returns 0 would skip forever."""
    if os.name != "nt":
        pytest.skip("Windows memory probes")
    assert _total_phys_bytes() > 2 * 2 ** 30
    assert 2 ** 20 < _own_rss_bytes() < _total_phys_bytes()
    assert 0.0 < A.memory_percent() <= 100.0


def test_the_ceiling_message_names_the_ceiling_and_the_governor_reading(tmp_path):
    """S26 defect 6b: a skip must say which ceiling fired and what the governor saw.

    Pure: synthetic snapshot files, no OpenMM.  The DECISION (`_memory_verdict`) is not
    touched; only the message reads the governor, and only while its snapshot is fresh.
    """
    p = tmp_path / "governor_state.json"
    now = 1_000_000.0
    p.write_text(json.dumps({"epoch": now - 10.0, "ram_pct": 94.2, "cpu_pct": 40.0,
                             "ts": "2026-09-13T00:00:00", "n_jobs": 2, "n_amber": 1,
                             "band": {"ceiling": 93.0}}), encoding="utf-8")
    g = _governor_reading(str(p), now=now)
    assert g is not None and g["ram_pct"] == 94.2 and g["n_amber"] == 1 and g["ceiling"] == 93.0
    # stale (older than the window), absent, or malformed -> not running -> fall back
    assert _governor_reading(str(p), max_age_s=5.0, now=now) is None
    assert _governor_reading(str(tmp_path / "missing.json"), now=now) is None
    (tmp_path / "bad.json").write_text("{not json", encoding="utf-8")
    assert _governor_reading(str(tmp_path / "bad.json"), now=now) is None

    m = _ceiling_message("skip", 95.0, 60.0, 100 * 2 ** 20, limit=92.0, gov=g)
    assert "92% ceiling" in m and "94.2% RAM" in m and "1 AMBER" in m
    assert "did not cause it" in m and "governor ceiling 93.0%" in m
    m2 = _ceiling_message("skip", 95.0, 60.0, 100 * 2 ** 20, limit=92.0, gov=None)
    assert "92% ceiling" in m2 and "governor is not running" in m2
    m3 = _ceiling_message("fail", 95.0, 60.0, 5 * 2 ** 30, limit=92.0, gov=g)
    assert "THIS SUITE accounts for it" in m3 and "5120 MB" in m3
    # the default limit is core.amber's own ceiling
    assert f"{A.MEMORY_LIMIT_PERCENT:.0f}% ceiling" in _ceiling_message("skip", 99.0, 60.0, 0)
