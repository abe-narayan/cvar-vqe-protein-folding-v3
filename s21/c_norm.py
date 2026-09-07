"""s21/c_norm.py -- BLOCK N: MEASURE BOTH COMPONENTS BEFORE ANYTHING IS COMBINED, AND THE GATES.

`s21/PREREG_C.md` section 1 declares the normalisation BEFORE any RMSD is read, and section 7
declares the gates.  This module produces both, and it is the only place in this lane where a
normalisation constant is manufactured.

WHAT IS MEASURED, per target, on the shipped top-75 REAL rebuilds (the domain the pipeline
actually sees) and separately at the c_land start set:

    numerical range and tails   mean median sd MAD IQR min max p1 p99 skew kurt frac>1e4 frac>1e6
    gradient scale              ||grad E|| at the medoid start, and its spread over 5 starts
    Hessian scale               ||H||_2 (WITHIN a potential -- has units) and the scale-free vector

for EACH potential SEPARATELY (`BRIEF` section 7: Legacy and AMBER separably evaluable always).

THE OBJECT BEING MEASURED, said plainly.  AMBER here is the **bare single point** --
`ConstrainedBox.energy_point`, no minimisation -- because a landscape needs a function of theta,
and `E o Relax_50` is not one (the relaxation leaves the torsion manifold).  The deployed
`H_AMBER` is `E o Relax_50` and that is a DIFFERENT operator; every table in this lane says which.

THE FOUR NORMALISATIONS.  Declared in `PREREG_C.md` section 1; `Nt` (asinh) is PRIMARY, and
`raw`, `Nz`, `Ng` are computed alongside so the choice is auditable.  `Nr` (rank->normal) is
declared REJECTED for continuation with its reason -- zero gradient a.e., undefined off-pool --
and is not implemented here for landscape use.

OPERATOR FORKS (`BRIEF` section 7 rule 0), with the alternative NOT taken named:
    FUNCTIONAL    TAKEN bare AMBER single point / NOT TAKEN `E o Relax_50` (not a function of
                  theta) and `E o Relax_1` (+inf on 42% of the register).
    BASIS         TAKEN the shipped top-75 real rebuilds / NOT TAKEN the k=8 lattice register.
    READOUT       none -- this module reads no RMSD at all, which is the point of running it
                  first.
    NORMALISATION TAKEN `Nt` PRIMARY / NOT TAKEN `raw`, `Nz`, `Ng` (all computed alongside) and
                  `Nr` (declared and rejected in advance: zero gradient a.e., undefined off-pool).
    NULL          none required -- no comparison is made here, only measurement.

    python -m s21.c_norm --gate      # GC21a/c/d
    python -m s21.c_norm --smoke     # 3 targets
    python -m s21.c_norm             # the 30-target c_land subset
    python -m s21.c_norm --report
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
from s15 import seed as SD                            # noqa: E402
from s18 import phys_lib as PL                        # noqa: E402
from s20 import c_land as CL                          # noqa: E402

#: the Sprint-20 subset, reused verbatim so every cross-sprint comparison is matched.
SUBSET = json.load(open(os.path.join(ROOT, "s20", "results", "c_land.json")))["config"]["subset"]
N_STARTS = 5


# ==========================================================================================
# THE NORMALISATION FAMILY.  Declared in PREREG_C.md section 1 before any RMSD was read.
# ==========================================================================================
class Norm:
    """An affine-or-asinh reparameterisation of one potential, with an EXACT derivative.

    Every member is C^1 and STRICTLY INCREASING, so by `PREREG_C.md` section 4 it leaves argmin,
    every CVaR tail SET and every critical point invariant.  That is a derivation, not a result,
    and it is why these transforms are safe to apply to a single Hamiltonian.
    """

    def __init__(self, kind, med=0.0, scale=1.0, gnorm=1.0):
        assert kind in ("raw", "Nz", "Ng", "Nt")
        self.kind = kind
        self.med = float(med)
        self.scale = float(max(scale, 1e-12))
        self.gnorm = float(max(gnorm, 1e-12))

    def __call__(self, E):
        E = np.asarray(E, float)
        if self.kind == "raw":
            return E
        if self.kind == "Nz":
            return (E - self.med) / self.scale
        if self.kind == "Ng":
            return (E - self.med) / self.gnorm
        z = (E - self.med) / self.scale
        return np.arcsinh(z)          # s*asinh(z)/s -- dimensionless, unit slope at z=0

    def dfdE(self, E):
        """f'(E) > 0 everywhere.  Used to map a raw gradient into the transformed one exactly."""
        E = np.asarray(E, float)
        if self.kind == "raw":
            return np.ones_like(E)
        if self.kind == "Nz":
            return np.full_like(E, 1.0 / self.scale)
        if self.kind == "Ng":
            return np.full_like(E, 1.0 / self.gnorm)
        z = (E - self.med) / self.scale
        return 1.0 / (self.scale * np.sqrt(1.0 + z * z))

    def d2fdE2(self, E):
        """f''(E).  Zero for the affine members; negative in the upper tail for `Nt`, which is
        the rank-one term of `PREREG_C.md` section 4.3."""
        E = np.asarray(E, float)
        if self.kind in ("raw", "Nz", "Ng"):
            return np.zeros_like(E)
        z = (E - self.med) / self.scale
        return -z / (self.scale ** 2 * (1.0 + z * z) ** 1.5)


def norms_for(stats_L, stats_A):
    """Build the four normalisations for one target from its own pool statistics.  No native
    information, no RMSD, no fitted constant beyond a robust pool centre and scale."""
    out = {}
    for kind in ("raw", "Nz", "Ng", "Nt"):
        out[kind] = {
            "legacy": Norm(kind, stats_L["median"], stats_L["mad_n"], stats_L["gnorm"]),
            "amber": Norm(kind, stats_A["median"], stats_A["mad_n"], stats_A["gnorm"]),
        }
    return out


# ==========================================================================================
def _desc(x):
    x = np.asarray(x, float)
    fin = np.isfinite(x)
    y = x[fin]
    med = float(np.median(y)) if len(y) else float("nan")
    mad = float(np.median(np.abs(y - med))) if len(y) else float("nan")
    q1, q3 = (np.percentile(y, [25, 75]) if len(y) else (np.nan, np.nan))
    sd = float(y.std(ddof=1)) if len(y) > 1 else float("nan")
    m3 = float(np.mean(((y - y.mean()) / sd) ** 3)) if len(y) > 2 and sd > 0 else float("nan")
    m4 = float(np.mean(((y - y.mean()) / sd) ** 4)) if len(y) > 3 and sd > 0 else float("nan")
    return {
        "n": int(len(x)), "n_finite": int(fin.sum()),
        "mean": float(y.mean()) if len(y) else float("nan"),
        "median": med, "sd": sd,
        "mad_n": float(1.4826 * mad),
        "iqr": float(q3 - q1),
        "min": float(y.min()) if len(y) else float("nan"),
        "max": float(y.max()) if len(y) else float("nan"),
        "p1": float(np.percentile(y, 1)) if len(y) else float("nan"),
        "p99": float(np.percentile(y, 99)) if len(y) else float("nan"),
        "skew": m3, "kurt": m4,
        "frac_gt_1e4": float((y > 1e4).mean()) if len(y) else float("nan"),
        "frac_gt_1e6": float((y > 1e6).mean()) if len(y) else float("nan"),
    }


#: SHARED-BOX GUARD FIRINGS.  `core.amber.memory_guard` refuses to open an OpenMM context above
#: 92% physical memory, because an OOM here kills a sibling workstream's run rather than only
#: this one.  It FIRED during Block N at target 16/30 with the box at 93% (an untouchable browser,
#: not a research process).  Retrying with backoff is legitimate -- the condition is transient and
#: external -- but the number of firings is a property of the run and is REPORTED, never hidden.
MEM_FIRINGS = {"n": 0, "waited_s": 0.0}


def with_mem_retry(fn, *a, tries=40, wait=45.0, **kw):
    import time as _t
    for k in range(tries):
        try:
            return fn(*a, **kw)
        except MemoryError as e:
            MEM_FIRINGS["n"] += 1
            MEM_FIRINGS["waited_s"] += wait
            print(f"    [memory_guard FIRED #{MEM_FIRINGS['n']}] {e}; waiting {wait:.0f}s",
                  flush=True)
            _t.sleep(wait)
    raise MemoryError(f"memory_guard fired {tries} times in a row; giving up")


def starts_for(pdb, PHI, PSI, W):
    """The IDENTICAL start rule s20/c_land.py used, reproduced exactly."""
    Pm = I.pairwise_rmsd(W)
    med = int(I.medoid(Pm))
    rng0 = SD.stable_rng(pdb, "s20C_land")
    others = [int(x) for x in rng0.permutation(len(W))[:N_STARTS] if int(x) != med]
    order = [med] + others[:N_STARTS - 1]
    return [np.concatenate([PHI[b], PSI[b]]) for b in order], order


def target_row(t):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    P = CL.Pot(seq)
    starts, order = starts_for(pdb, PHI, PSI, W)
    active = np.arange(1, 2 * n)
    TH = np.column_stack([PHI, PSI])

    e_l = np.asarray(P.legacy(TH), float)
    e_a = np.asarray(P.amber(TH), float)
    sL, sA = _desc(e_l), _desc(e_a)

    #: gradient scale at every start (NOT only the medoid -- the spread is the point)
    gl, ga = [], []
    for th in starts:
        gl.append(float(np.linalg.norm(P.grad_fd("legacy", th)[active])))
        ga.append(float(np.linalg.norm(P.amber_grad_chain(th)[1][active])))
    sL["gnorm"] = float(gl[0]); sA["gnorm"] = float(ga[0])
    sL["gnorm_starts"] = gl;    sA["gnorm_starts"] = ga
    sL["gnorm_cv"] = float(np.std(gl) / max(np.mean(gl), 1e-30))
    sA["gnorm_cv"] = float(np.std(ga) / max(np.mean(ga), 1e-30))

    #: Hessian scale at the medoid start, WITHIN each potential (units) plus the scale-free vector
    for which, s in (("legacy", sL), ("amber", sA)):
        Hm, e0 = P.hess_fd(which, starts[0], active=active)
        m = CL.spectrum_metrics(Hm)
        s["E0"] = float(e0)
        s["hess_norm2"] = float(m["lam_max_abs"])
        s["spectrum"] = {k: m[k] for k in CL.METRICS}

    #: THE RAW-UNITS CROSSOVER lambda* -- section 0.3 of the pre-registration.  At which lambda
    #: does the AMBER term first dominate the mixed gradient, in RAW units?
    lam_star = [float(a / (a + b)) for a, b in zip(gl, ga)]

    P.close()
    return {"pdb": pdb, "n": n, "fold": int(t["fold"]),
            "legacy": sL, "amber": sA,
            "lam_star_raw": lam_star,
            "lam_star_raw_med": float(np.median(lam_star)),
            "pool_m": int(len(TH))}


# ==========================================================================================
# GATES
# ==========================================================================================
def gate(n_targets=3, verbose=True):
    """GC21a (instrument reproduces Sprint 20), GC21c (Nt monotone, MAD positive), GC21d (seal)."""
    out = {"GC21a": [], "GC21c": {}, "GC21d": {}, "firings": {}}
    old = {r["pdb"]: r for r in
           json.load(open(os.path.join(ROOT, "s20", "results", "c_land.json")))["rows"]}
    tg = {t["pdb"]: t for t in I.targets()}
    fired_a = 0
    for pdb in SUBSET[:n_targets]:
        t = tg[pdb]; seq = t["seq"]; n = int(t["n"])
        W, PHI, PSI, u = top75_windows(pdb)
        W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
        P = CL.Pot(seq)
        starts, _ = starts_for(pdb, PHI, PSI, W)
        active = np.arange(1, 2 * n)
        r = {"pdb": pdb}
        for which in ("legacy", "amber"):
            Hm, e0 = P.hess_fd(which, starts[0], active=active)
            m = CL.spectrum_metrics(Hm)
            ref = old[pdb]["start"][which]
            d = {k: float(abs(m[k] - ref[k])) for k in CL.METRICS}
            d["E0"] = float(abs(e0 - ref["E0"]))
            r[which] = d
            if max(d.values()) > 1e-9:
                fired_a += 1
        P.close()
        out["GC21a"].append(r)
        if verbose:
            print(f"GC21a {pdb}: max|delta| legacy {max(r['legacy'].values()):.3e} "
                  f"amber {max(r['amber'].values()):.3e}", flush=True)

    #: GC21c -- Nt strict monotonicity, exactly, on a wide synthetic range plus the real one.
    fired_c = 0
    x = np.sort(np.concatenate([
        np.linspace(-2000, 5000, 4001),
        np.geomspace(1e3, 1e24, 4000),
        -np.geomspace(1e0, 1e6, 1000)[::-1]]))
    nt = Norm("Nt", med=16062.0, scale=1e4)
    y = nt(x)
    mono = bool(np.all(np.diff(y) > 0))
    dpos = bool(np.all(nt.dfdE(x) > 0))
    if not (mono and dpos):
        fired_c += 1
    out["GC21c"] = {"monotone_on_synthetic_range": mono, "derivative_positive": dpos,
                    "n_pairs_checked": int(len(x) - 1),
                    "range": [float(x.min()), float(x.max())],
                    "max_finite_output": float(np.max(np.abs(y)))}

    #: GC21d -- the benchmark seal, as an integrity certificate that needs NO read.
    man = os.path.join(ROOT, "results", "benchmark_manifest.json")
    h = hashlib.sha256()
    if os.path.exists(man):
        with open(man, "rb") as fh:
            for blk in iter(lambda: fh.read(1 << 20), b""):
                h.update(blk)
        out["GC21d"] = {"path": "results/benchmark_manifest.json", "sha256": h.hexdigest(),
                        "bytes": os.path.getsize(man), "read_by_this_lane": False}
    else:
        out["GC21d"] = {"path": man, "exists": False}

    out["firings"] = {"GC21a_reproduction": int(fired_a), "GC21c_monotonicity": int(fired_c)}
    out["passed"] = bool(fired_a == 0 and fired_c == 0)
    if verbose:
        print(f"GC21c: Nt monotone over {out['GC21c']['n_pairs_checked']} pairs spanning "
              f"[{out['GC21c']['range'][0]:.1e}, {out['GC21c']['range'][1]:.1e}]  -> "
              f"{out['GC21c']['monotone_on_synthetic_range']}")
        print(f"GC21d: benchmark_manifest.json sha256 {out['GC21d'].get('sha256','(absent)')} "
              f"-- recorded WITHOUT reading its content")
        print(f"GATES {'PASS' if out['passed'] else 'FAIL'}  firings {out['firings']}")
    json.dump(out, open(os.path.join(RESULTS, "c_gate.json"), "w"), indent=1)
    return out


# ==========================================================================================
def run(subset=None, out="c_norm.json"):
    sub = SUBSET if subset is None else subset
    tg = {t["pdb"]: t for t in I.targets()}
    path = os.path.join(RESULTS, out)
    rows, done = [], set()
    if os.path.exists(path):
        prev = json.load(open(path))
        rows = [r for r in prev.get("rows", []) if r["pdb"] in sub]
        done = {r["pdb"] for r in rows}
    t0 = time.time()
    for pdb in sub:
        if pdb in done:
            continue
        rows.append(with_mem_retry(target_row, tg[pdb]))
        r = rows[-1]
        print(f"  {len(rows)}/{len(sub)} {pdb}  E_L med {r['legacy']['median']:+.2f}  "
              f"E_A med {r['amber']['median']:+.3e}  E_A max {r['amber']['max']:.3e}  "
              f"lam*_raw {r['lam_star_raw_med']:.3e}  ({time.time()-t0:.0f}s)", flush=True)
        _write(rows, sub, path)
    _write(rows, sub, path)
    report()
    return rows


def _write(rows, sub, path):
    cfg = {"subset": list(sub), "N_STARTS": N_STARTS, "H_GRAD": CL.H_GRAD, "H_HESS": CL.H_HESS,
           "amber_object": "bare single point (ConstrainedBox.energy_point), NO minimisation",
           "norms": ["raw", "Nz", "Ng", "Nt"], "primary": "Nt",
           "rejected_for_continuation": {"Nr": "zero gradient a.e.; undefined off-pool"}}
    cfg["memory_guard_firings"] = int(MEM_FIRINGS["n"])
    obj = {"rows": rows, "config": cfg, "memory_guard_firings": int(MEM_FIRINGS["n"]),
           "memory_guard_wait_s": float(MEM_FIRINGS["waited_s"]),
           "cfg_hash": hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16],
           "n_rows": len(rows), "n_expected": len(sub),
           "complete": len(rows) >= len(sub)}
    json.dump(obj, open(path, "w"))


def report(out="c_norm.json"):
    o = json.load(open(os.path.join(RESULTS, out)))
    rows = o["rows"]
    tag = "" if o["complete"] else f"   *** PARTIAL n={len(rows)} of {o['n_expected']} ***"
    L = "=" * 104
    print(L)
    print(f"BLOCK N -- COMPONENT MEASUREMENT BEFORE ANY COMBINATION   n = {len(rows)} targets{tag}")
    print("  domain: the shipped top-75 REAL rebuilds per target.  AMBER = BARE SINGLE POINT, no")
    print("  minimisation.  The deployed H_AMBER is E o Relax_50 and is a DIFFERENT operator.")
    print(L)
    g = lambda p, k: np.array([r[p][k] for r in rows], float)          # noqa: E731
    print(f"\n{'statistic':<16}{'LEGACY (arb. units)':>26}{'AMBER (kcal/mol)':>26}")
    for k, lab in (("median", "pool median"), ("mad_n", "pool 1.4826*MAD"),
                   ("sd", "pool sd"), ("iqr", "pool IQR"),
                   ("min", "pool min"), ("max", "pool max"),
                   ("p1", "pool p1"), ("p99", "pool p99"),
                   ("skew", "skew"), ("kurt", "kurtosis"),
                   ("frac_gt_1e4", "frac > 1e4"), ("frac_gt_1e6", "frac > 1e6"),
                   ("gnorm", "||grad|| at medoid"), ("gnorm_cv", "||grad|| CV over starts"),
                   ("hess_norm2", "||Hess||_2"), ("E0", "E at medoid start")):
        a, b = g("legacy", k), g("amber", k)
        print(f"{lab:<16}{np.median(a):>26.4g}{np.median(b):>26.4g}")
    print("\n  (medians over targets.  Raw scales are NOT comparable across the two columns --")
    print("   that is the entire reason this block exists.)")

    ls = g("legacy", "median") * 0 + np.array([r["lam_star_raw_med"] for r in rows], float)
    print(f"\n  RAW-UNIT CROSSOVER  lambda* = ||grad E_L|| / (||grad E_L|| + ||grad E_A||)")
    print(f"    median {np.median(ls):.3e}   min {ls.min():.3e}   max {ls.max():.3e}   "
          f"spread {ls.max()/max(ls.min(),1e-30):.1f}x")
    print(f"    -> in RAW units a uniform lambda grid is AMBER-dominated above lambda ~ "
          f"{np.median(ls):.2e}.")

    print(f"\n  SCALE-FREE SPECTRUM at the medoid start (reproducing s20 L12):")
    print(f"    {'metric':<16}{'legacy':>12}{'amber':>12}")
    for k in CL.METRICS:
        a = np.array([r["legacy"]["spectrum"][k] for r in rows], float)
        b = np.array([r["amber"]["spectrum"][k] for r in rows], float)
        print(f"    {k:<16}{np.median(a):>12.4g}{np.median(b):>12.4g}")

    #: what each normalisation DOES to the tail -- the audit the pre-registration promised
    print(f"\n  WHAT EACH DECLARED NORMALISATION DOES TO THE WORST REAL AMBER VALUE IN THIS SET:")
    imax = int(np.argmax(g("amber", "max")))
    r = rows[imax]
    med, sc, gn = r["amber"]["median"], r["amber"]["mad_n"], r["amber"]["gnorm"]
    worst = r["amber"]["max"]
    for kind in ("raw", "Nz", "Ng", "Nt"):
        nz = Norm(kind, med, sc, gn)
        print(f"    {kind:<5}{float(nz(np.array([worst]))[0]):>18.4g}   "
              f"(pool median maps to {float(nz(np.array([med]))[0]):+.4g})")
    print(f"    target {r['pdb']}: raw worst {worst:.4g}, pool median {med:.4g}, MAD_n {sc:.4g}")
    print(L)


if __name__ == "__main__":
    a = sys.argv[1:]
    if "--gate" in a:
        gate()
    elif "--report" in a:
        report()
    elif "--smoke" in a:
        run(subset=SUBSET[:3], out="_SMOKE_c_norm.json")
    else:
        run()
