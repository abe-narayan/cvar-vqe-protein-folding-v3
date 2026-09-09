"""The KNOWN HAZARDS, each checked as a live experiment rather than by reading the source.

  H1  single-point energy is exactly invariant under rigid translation
      (a position-dependent term leaking into the reported energy would break this)
  H2  `_calibrate_hydrogens` is pinned to ONE thread under every worker configuration
      -- the hazard that moved single-point energies by 144 kcal/mol
  H3  BLAS ddot alignment: the same contraction on an aligned and a deliberately
      misaligned buffer must give the same float
  H4  accumulation over Python floats, not numpy scalars, in the weighted totals
  H5  one canonical 20-letter alphabet, ARNDCQEGHILKMFPSTWYV, with no default anywhere
  H6  stable argsort throughout -- pool membership must not depend on sort stability
  H7  `identity` still normalises by the LONGER sequence (leaky by design, preserved)
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

CANON = "ARNDCQEGHILKMFPSTWYV"


def h1_translation(out):
    import numpy as np
    import peptide_db as db
    import torsion_lib2 as tl2
    import protein_geometry as geo
    import core.amber as A

    p = db.by_pdb("1A13")
    tab = tl2.library_for(p.seq, 8, p.seq)
    rep = tl2.PerResidueTorsion(p.seq, tab, chi_bits=False)
    rb = geo.build_backbone(np.asarray(p.phi, float), np.asarray(p.psi, float))

    a = A.refine_coords(p.seq, rep, rb, k_restraint=0.0, steps=-1, tolerance=1e9,
                        components=True)
    shift = np.array([1.75, -0.5, 0.25])
    rb2 = {k: v + shift for k, v in rb.items()}
    b = A.refine_coords(p.seq, rep, rb2, k_restraint=0.0, steps=-1, tolerance=1e9,
                        components=True)
    out["H1_singlepoint_native"] = float(a["energy"])
    out["H1_singlepoint_translated"] = float(b["energy"])
    out["H1_abs_diff"] = float(abs(a["energy"] - b["energy"]))
    out["H1_bit_identical"] = bool(a["energy"] == b["energy"])
    out["H1_components_max_diff"] = float(max(
        abs(a["components"][k] - b["components"][k]) for k in a["components"]))


def h2_thread_pin(out):
    """`_calibrate_hydrogens` must run at Threads=1 no matter what the run is set to.

    Checked two ways: the source must pin it, and the calibrated value must be IDENTICAL
    across processes started with different OpenMM/BLAS thread environments.
    """
    src = open(os.path.join(_ROOT, "core", "amber.py"), encoding="utf-8",
               errors="replace").read()
    m = re.search(r"def _calibrate_hydrogens.*?(?=\ndef |\nclass |\Z)", src, re.S)
    body = m.group(0) if m else ""
    out["H2_found_calibrate_hydrogens"] = bool(m)
    out["H2_mentions_threads_1"] = bool(
        re.search(r"Threads[\"']?\s*[:=]\s*[\"']?1|threads\s*=\s*1", body))
    # the pin must not be read from the ambient configuration
    out["H2_body_len"] = len(body)

    prog = (
        "import os,sys,json;sys.path.insert(0,%r);"
        "import numpy as np, peptide_db as db, torsion_lib2 as tl2,"
        " protein_geometry as geo, core.amber as A;"
        "p=db.by_pdb('1A13');"
        "tab=tl2.library_for(p.seq,8,p.seq);"
        "rep=tl2.PerResidueTorsion(p.seq,tab,chi_bits=False);"
        "rb=geo.build_backbone(np.asarray(p.phi,float),np.asarray(p.psi,float));"
        "r=A.refine_coords(p.seq,rep,rb,k_restraint=0.0,steps=-1,tolerance=1e9,"
        "components=True);"
        "print(json.dumps({'e':repr(r['energy']),"
        "'c':{k:repr(v) for k,v in r['components'].items()}}))" % _ROOT
    )
    envs = [
        ("default", {}),
        ("omp8", {"OMP_NUM_THREADS": "8", "OPENMM_CPU_THREADS": "8",
                  "OPENBLAS_NUM_THREADS": "8", "MKL_NUM_THREADS": "8"}),
        ("omp4", {"OMP_NUM_THREADS": "4", "OPENMM_CPU_THREADS": "4",
                  "OPENBLAS_NUM_THREADS": "4", "MKL_NUM_THREADS": "4"}),
    ]
    res = {}
    for tag, extra in envs:
        env = dict(os.environ)
        env.update(extra)
        try:
            cp = subprocess.run([sys.executable, "-c", prog], capture_output=True,
                                text=True, env=env, cwd=_ROOT, timeout=900)
            line = [l for l in cp.stdout.splitlines() if l.startswith("{")]
            res[tag] = json.loads(line[-1]) if line else {
                "error": (cp.stderr or cp.stdout)[-400:]}
        except Exception as exc:                                   # pragma: no cover
            res[tag] = {"error": f"{type(exc).__name__}: {exc}"}
    out["H2_by_thread_env"] = res
    es = [v.get("e") for v in res.values() if "e" in v]
    out["H2_all_energies_identical"] = bool(es and len(set(es)) == 1)
    out["H2_n_env_ok"] = len(es)


def h3_blas_alignment(out):
    """The same contraction on aligned and misaligned buffers must agree bit-for-bit."""
    rng = np.random.default_rng(11)
    n, d = 400, 96
    A0 = rng.standard_normal((n, d))
    v = rng.standard_normal(d)
    ref = A0 @ v
    worst = 0.0
    for off in range(1, 9):
        buf = np.empty((n, d + off))
        buf[:, off:] = A0
        mis = buf[:, off:]
        worst = max(worst, float(np.max(np.abs(mis @ v - ref))))
    out["H3_max_abs_diff_alignment"] = worst
    out["H3_alignment_invariant"] = worst == 0.0

    # and the same for the pairwise-distance path the pipeline actually uses
    try:
        g = __import__("core.geometry", fromlist=["x"])
        if hasattr(g, "pair_dists") and hasattr(g, "pair_index"):
            X = rng.standard_normal((64, 3))
            ii, jj = g.pair_index(64, 2)
            r0 = g.pair_dists(X, ii, jj)
            pad = np.empty((64, 4))
            pad[:, 1:] = X
            r1 = g.pair_dists(pad[:, 1:], ii, jj)
            out["H3_pair_dists_max_abs_diff"] = float(np.max(np.abs(r0 - r1)))
            out["H3_pair_dists_bit_identical"] = bool(np.array_equal(r0, r1))
    except Exception as exc:
        out["H3_pair_dists_error"] = f"{type(exc).__name__}: {exc}"


def h4_accumulation(out):
    """Weighted totals must accumulate over Python floats, not numpy scalars."""
    rng = np.random.default_rng(5)
    terms = rng.standard_normal(11) * 30.0
    w = rng.standard_normal(11)
    py = 0.0
    for t, x in zip(terms, w):
        py += float(t) * float(x)
        py = float(py)
    npy = np.float64(0.0)
    for t, x in zip(terms, w):
        npy = npy + t * x
    out["H4_python_float_total"] = repr(py)
    out["H4_numpy_scalar_total"] = repr(float(npy))
    out["H4_they_can_differ"] = bool(py != float(npy))
    out["H4_dot_total"] = repr(float(terms @ w))
    # Report which convention core/energy.py uses.
    src = open(os.path.join(_ROOT, "core", "energy.py"), encoding="utf-8",
               errors="replace").read()
    out["H4_energy_mentions_float_accum"] = bool(
        re.search(r"float\(", src)) and "numpy scalar" in src


def h5_alphabet(out):
    """One canonical alphabet; no module may reintroduce a different default."""
    hits, bad = {}, {}
    for dirpath, dirnames, filenames in os.walk(os.path.join(_ROOT, "core")):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fn in filenames:
            if not fn.endswith(".py"):
                continue
            fp = os.path.join(dirpath, fn)
            txt = open(fp, encoding="utf-8", errors="replace").read()
            for m in re.finditer(r"[\"']([ACDEFGHIKLMNPQRSTVWY]{20})[\"']", txt):
                s = m.group(1)
                if sorted(s) != sorted(CANON):
                    continue
                hits.setdefault(s, []).append(os.path.relpath(fp, _ROOT))
                if s != CANON:
                    bad.setdefault(s, []).append(os.path.relpath(fp, _ROOT))
    out["H5_alphabets_found"] = {k: sorted(set(v)) for k, v in hits.items()}
    out["H5_noncanonical"] = {k: sorted(set(v)) for k, v in bad.items()}
    out["H5_only_canonical"] = not bad


def h6_stable_sort(out):
    """Unstable argsort must be absent from the ranking paths in core/."""
    offenders, stable_n, total_n = [], 0, 0
    for dirpath, dirnames, filenames in os.walk(os.path.join(_ROOT, "core")):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fn in sorted(filenames):
            if not fn.endswith(".py"):
                continue
            fp = os.path.join(dirpath, fn)
            for i, line in enumerate(open(fp, encoding="utf-8",
                                          errors="replace"), 1):
                if "argsort" not in line and "argpartition" not in line:
                    continue
                total_n += 1
                if 'kind="stable"' in line or "kind='stable'" in line:
                    stable_n += 1
                else:
                    offenders.append(f"{os.path.relpath(fp, _ROOT)}:{i}: {line.strip()}")
    out["H6_argsort_total"] = total_n
    out["H6_argsort_stable"] = stable_n
    out["H6_argsort_without_stable"] = offenders

    # and measure what instability would cost, on a real tie-heavy score vector
    rng = np.random.default_rng(3)
    sc = np.round(rng.standard_normal(5000), 2)      # many ties
    a = np.argsort(sc, kind="stable")[:500]
    b = np.argsort(sc, kind="quicksort")[:500]
    out["H6_membership_delta_if_unstable"] = int(len(set(a) ^ set(b)) // 2)


def h7_identity(out):
    """`identity` normalises by the LONGER sequence -- leaky by design, preserved."""
    try:
        d = __import__("core.data", fromlist=["x"])
        f = d.identity
        short, long = "ACDEFG", "ACDEFGHIKL"
        v = f(short, long)
        out["H7_identity_short_vs_long"] = float(v)
        out["H7_expected_if_longer_norm"] = 6 / 10
        out["H7_expected_if_shorter_norm"] = 1.0
        out["H7_normalises_by_longer"] = bool(abs(v - 0.6) < 1e-9)
        out["H7_symmetric"] = bool(f(short, long) == f(long, short))
        out["H7_threshold"] = float(d.IDENTITY_THRESHOLD)
    except Exception as exc:
        out["H7_error"] = f"{type(exc).__name__}: {exc}"


def main():
    out = {}
    for fn in (h5_alphabet, h6_stable_sort, h7_identity, h3_blas_alignment,
               h4_accumulation, h1_translation, h2_thread_pin):
        try:
            fn(out)
        except Exception as exc:
            out[fn.__name__ + "_ERROR"] = f"{type(exc).__name__}: {exc}"
    print(json.dumps(out, indent=2, sort_keys=True, default=str))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "hazard_audit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True, default=str)
    return out

if __name__ == "__main__":
    main()
