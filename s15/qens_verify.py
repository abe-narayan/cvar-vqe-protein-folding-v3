"""SPRINT 15 / QENS -- verification before any claim.  Five independent checks.

V1  the instrument reproduces its five pinned constants.
V2  `GEnum` + `Struct` reproduce `s15.qgeom_cvar.Struct` EXACTLY on the nine n=9 targets
    (identical `full` index vector, identical ORACLE rmsd column), and the precomputed CA
    table reproduces the stored ORACLE rmsd column to machine precision on all 19.
V3  the sampled-sd objective of `qens_lib.hamil_sub` is RANK-IDENTICAL to the cached
    full-tabulation objective `s14.vqe_hamil.combine(*tabulate(pdb), 0.25)[full]` on all nine
    n=9 targets.  If it is, every VQE run in this workstream is bit-identical to what
    `qgeom_ens` would have produced on those targets, and the 19-target extension is on one
    code path.
V4  the CVaR-VQE run reproduces bit-identically from the same seed (determinism).
V5  timing, so the experiment can be sized before it is launched.

    python -m s15.qens_verify
"""
from __future__ import annotations

import time

import numpy as np

from s15 import qens_lib as Q
from s15 import qgeom_lib as G


def v1_instrument():
    """Run `python -m s12.instrument` as a subprocess and check its printed constants.

    NOTE, from the Phase 0 audit: this is a CACHE READ, not a reproduction of the pipeline.
    It is run here to confirm the shared instrument is the pinned one, nothing more.
    """
    import re
    import subprocess
    import sys
    print("=" * 100)
    print("V1. INSTRUMENT (cache read, per the Phase 0 audit -- not a pipeline reproduction)")
    print("=" * 100)
    exp = {"shipped": 3.4540004952559396, "pool_best": 1.7108244199364904,
           "top75_best": 2.3061526409453816, "synthesis_fit": 3.2040761603809194}
    p = subprocess.run([sys.executable, "-m", "s12.instrument"], cwd=Q.ROOT,
                       capture_output=True, text=True, timeout=900)
    txt = p.stdout + p.stderr
    got, ok = {}, True
    for k, v in exp.items():
        m = re.search(rf"{k}\D{{0,12}}([0-9]+\.[0-9]+)", txt)
        g = float(m.group(1)) if m else None
        got[k] = g
        hit = g is not None and abs(g - v) < 1e-9
        ok &= hit
        print(f"  {k:16s} expected {v:.10f}  got {g if g is None else f'{g:.10f}'}  "
              f"{'OK' if hit else 'MISMATCH'}")
    if not ok:
        print("  raw tail:", txt.strip().splitlines()[-6:])
    Q.ck("verify", "V1_instrument", {"ok": bool(ok), "got": got})
    return {"ok": bool(ok), "got": got}


def v2_struct():
    print()
    print("=" * 100)
    print("V2. STRUCT / CA TABLE")
    print("=" * 100)
    from s15.qgeom_cvar import Struct as OldStruct
    out = {}
    for pdb in Q.TARGETS9:
        a = Q.Struct(pdb)
        b = OldStruct(pdb, Q.SUB)
        same_full = bool(np.array_equal(a.full, b.full))
        same_rmsd = float(np.abs(a.rmsd - b.rmsd).max())
        out[pdb] = {"full_identical": same_full, "rmsd_maxdiff": same_rmsd,
                    "ca_rebuild_maxerr": a.rebuild_maxerr, "n": a.e.n}
        print(f"  {pdb}  full identical {same_full}  rmsd maxdiff {same_rmsd:.2e}  "
              f"CA rebuild maxerr {a.rebuild_maxerr:.2e}")
    for pdb in Q.TARGETS10:
        a = Q.Struct(pdb)
        out[pdb] = {"full_identical": None, "rmsd_maxdiff": None,
                    "ca_rebuild_maxerr": a.rebuild_maxerr, "n": a.e.n}
        print(f"  {pdb}  n={a.e.n}  M={a.M}  CA rebuild maxerr {a.rebuild_maxerr:.2e}  "
              f"ORACLE rmsd range {a.rmsd.min():.3f}-{a.rmsd.max():.3f}")
    Q.ck("verify", "V2_struct", out)
    return out


def v3_objective():
    print()
    print("=" * 100)
    print("V3. qens_lib.hamil_sub vs THE ORIGINAL qgeom_ens LINE (nine n=9 targets)")
    print("=" * 100)
    from s14.vqe_hamil import tabulate, combine
    out = {}
    for pdb in Q.TARGETS9:
        st = Q.Struct(pdb)
        pr, ds = tabulate(pdb)
        E_ref = G.V.uniformise(combine(pr, ds, 0.25)[st.full])
        E_new = Q.hamil_sub(st, w=0.25)
        ident = bool(np.array_equal(np.argsort(E_ref, kind="mergesort"),
                                    np.argsort(E_new, kind="mergesort")))
        out[pdb] = {"rank_identical": ident,
                    "max_abs_diff": float(np.abs(E_ref - E_new).max()),
                    "spearman": float(G.V.spearman(E_ref, E_new))}
        print(f"  {pdb}  rank-identical {ident}   max|dE| {out[pdb]['max_abs_diff']:.3e}   "
              f"rho {out[pdb]['spearman']:.9f}")
    Q.ck("verify", "V3_objective", out)
    return out


def v4_determinism():
    print()
    print("=" * 100)
    print("V4. DETERMINISM of the CVaR-VQE run")
    print("=" * 100)
    from s15.qgeom_cvar import make_circ, run_vqe
    st = Q.Struct("1CS9")
    E = Q.hamil_sub(st)
    circ = make_circ(st.n_qubits)
    a = run_vqe(circ, E, 0.05, 60, seed=3, lowmem=True)
    b = run_vqe(circ, E, 0.05, 60, seed=3, lowmem=True)
    d = float(np.abs(a["p_final"] - b["p_final"]).max())
    print(f"  max |p_final - p_final'| over a repeated run: {d:.3e}")
    Q.ck("verify", "V4_determinism", {"max_p_diff": d})
    return d


def v5_timing():
    print()
    print("=" * 100)
    print("V5. TIMING")
    print("=" * 100)
    from s15.qgeom_cvar import make_circ, run_vqe
    t0 = time.time()
    st = Q.Struct("2MJQ")
    t_struct = time.time() - t0
    t0 = time.time()
    E = Q.hamil_sub(st, cache=False)
    t_obj = time.time() - t0
    circ = make_circ(st.n_qubits)
    t0 = time.time()
    r = run_vqe(circ, E, 0.05, 200, seed=0, lowmem=True)
    t_vqe = time.time() - t0
    rng = np.random.default_rng(0)
    idx = rng.choice(len(E), size=2048, p=r["p_final"])
    t0 = time.time()
    ro = Q.readouts(st, idx, np.random.default_rng(0), dmatch=50)
    t_ro = time.time() - t0
    t0 = time.time()
    ai, used = Q.anneal_ensemble(E, len(st.res), st.e.k, 2048, np.random.default_rng(0))
    t_sa = time.time() - t0
    print(f"  Struct build (n=10)      {t_struct:7.2f}s")
    print(f"  objective (sampled sd)   {t_obj:7.2f}s   [cached thereafter]")
    print(f"  CVaR-VQE 200 iters       {t_vqe:7.2f}s")
    print(f"  one arm's readouts       {t_ro:7.2f}s")
    print(f"  SA 2048 evals            {t_sa:7.2f}s  (used {used})")
    print(f"  -> per cell with 6 arms  ~{t_vqe + 6 * t_ro + t_sa:7.2f}s")
    Q.ck("verify", "V5_timing", {"struct": t_struct, "objective": t_obj, "vqe200": t_vqe,
                                 "readouts": t_ro, "sa2048": t_sa})
    return t_vqe, t_ro


def main():
    G.wait_mem(0.8, "qens_verify")
    v1_instrument()
    v2_struct()
    v3_objective()
    v4_determinism()
    v5_timing()
    print("\nwritten -> s15/results/qens_verify.json")


if __name__ == "__main__":
    main()
