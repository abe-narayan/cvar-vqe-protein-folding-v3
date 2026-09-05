"""`core.bench` -- the measurement instrument for the speedup claim.

WHY THIS FILE EXISTS AND WHAT IT REFUSES TO DO
==============================================
A speedup claim without a committed baseline, measured on the same scientific experiment,
is worthless.  So this harness:

  * runs the COMPLETE pipeline at PRODUCTION settings on a DECLARED target manifest.  The
    headline is the 126-target tuning instrument.  Smaller manifests exist, are labelled
    `dev_mode`, and are refused as a headline;
  * emits per-stage wall clock -- retrieval, distogram scoring, filtering, coordinate
    averaging, projection, AMBER -- plus I/O, serialisation, process startup and import
    time, because a "10x on the hot loop" that leaves a 40-second import untouched is not
    a 10x;
  * records CPU%, RAM peak and mean, and the worker/thread configuration alongside EVERY
    timing, so a number that was really bought with oversubscription or with the box
    swapping is visible as such;
  * writes machine-readable results to `bench_results/*.json`.  Speedups are computed by
    `--compare` from two real runs on disk.  Nothing here estimates anything.

THE TWO ARMS
============
`--arm baseline` runs the `s9/final.py` path: its own `library_members`, its own
`windows_all`, its own full K x K pairwise matrix, its own `synthesise`, its own AMBER
block, one target at a time in one process, with the float32 `npz` round trip it puts
between its CLI stages actually performed.  It calls `s9.final`'s functions; it does not
reimplement them.

`--arm optimised` runs `core.pipeline`.

Both arms write per-target checkpoints and are resumable, and both report `n_resumed` --
a run that inherited work from disk is not a cold measurement and the results file says so.

    python -m core.bench --arm baseline  --manifest tuning126 --fresh
    python -m core.bench --arm optimised --manifest tuning126 --workers 6 --fresh
    python -m core.bench --compare bench_results/baseline_tuning126.json \\
                                   bench_results/optimised_tuning126.json
    python -m core.bench --sweep-workers 1,2,4,6,8 --manifest smoke24
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import threading
import time
from dataclasses import asdict, replace

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np                                                       # noqa: E402

import core                                                              # noqa: E402
from core import pipeline as P                                           # noqa: E402

RESULTS = P.RESULTS
BASE_CACHE = os.path.join(P.CACHE_ROOT, "baseline_s9")


# ============================================================ resource sampling
class Resources:
    """System CPU% and RAM, sampled on a background thread while a run is in flight.

    `GetSystemTimes` for CPU (idle/kernel/user deltas -- whole box, which is what "target
    90-95% CPU" means when the work is spread over worker processes) and
    `GlobalMemoryStatusEx` for RAM.  Deliberately NOT `Get-CimInstance`: on a loaded box
    that costs minutes per call, so the instrument would become the load.
    """

    def __init__(self, period=0.25):
        self.period = period
        self.cpu, self.ram_pct, self.ram_used = [], [], []
        self._stop = threading.Event()
        self._th = None
        self.total_ram = P.mem_bytes()[0]

    @staticmethod
    def _systimes():
        import ctypes
        from ctypes import wintypes

        class FT(ctypes.Structure):
            _fields_ = [("lo", wintypes.DWORD), ("hi", wintypes.DWORD)]

        def v(f):
            return (f.hi << 32) | f.lo
        i, k, u = FT(), FT(), FT()
        ctypes.windll.kernel32.GetSystemTimes(ctypes.byref(i), ctypes.byref(k),
                                              ctypes.byref(u))
        return v(i), v(k), v(u)

    def _loop(self):
        prev = self._systimes()
        while not self._stop.wait(self.period):
            cur = self._systimes()
            di = cur[0] - prev[0]
            dk = cur[1] - prev[1]
            du = cur[2] - prev[2]
            prev = cur
            busy = (dk + du) - di            # kernel time includes idle
            tot = dk + du
            if tot > 0:
                self.cpu.append(100.0 * busy / tot)
            total, avail, pct = P.mem_bytes()
            self.ram_pct.append(float(pct))
            self.ram_used.append(float(total - avail))

    def __enter__(self):
        self._th = threading.Thread(target=self._loop, daemon=True)
        self._th.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        if self._th:
            self._th.join(timeout=2.0)
        return False

    def report(self):
        def s(v, f=lambda x: x):
            if not v:
                return None
            a = np.asarray(v, float)
            return {"mean": float(f(a.mean())), "peak": float(f(a.max())),
                    "min": float(f(a.min()))}
        return {"n_samples": len(self.cpu),
                "sample_period_s": self.period,
                "cpu_pct": s(self.cpu),
                "ram_pct": s(self.ram_pct),
                "ram_used_gb": s(self.ram_used, lambda x: x / 2 ** 30),
                "total_ram_gb": self.total_ram / 2 ** 30}


# ============================================================ startup / import cost
#: The heavy closure a worker really touches once it starts doing work: torch (the
#: distogram predictor), OpenMM (the relaxation) and the peptide database.  Several of
#: these are DEFERRED imports inside the stage functions, so without this measurement they
#: would be silently charged to whichever stage happened to touch them first -- and a
#: "0.12 s import" would be a fiction.
_HEAVY = ("import numpy, torch, openmm; "
          "import peptide_db, protein_geometry, distogram, amber_refine, torsion_lib2; "
          "from s7 import audit, debias; from s8 import consensus2, project, inband")


def measure_startup(module="core.pipeline"):
    """Process startup and import time, measured with real subprocesses.

    `startup_s` is a bare interpreter.  `import_s` is the same interpreter plus the arm's
    own module, minus the bare cost.  `heavy_import_s` adds the deferred closure the
    stages pull in on first use.  All three are PER-PROCESS costs, so they are the tax a
    parallel run pays once per worker and a serial run pays once -- which is exactly why
    they are reported apart from the stage clocks and never folded into a speedup.
    """
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = env["MKL_NUM_THREADS"] = "1"

    def t(code, n=3):
        best = float("inf")
        for _ in range(n):
            t0 = time.perf_counter()
            subprocess.run([sys.executable, "-c", code], cwd=_ROOT, env=env,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           check=False)
            best = min(best, time.perf_counter() - t0)
        return best
    pre = f"import sys; sys.path.insert(0, {_ROOT!r}); "
    bare = t("pass")
    full = t(pre + f"import {module}")
    heavy = t(pre + f"import {module}; " + _HEAVY, n=2)
    return {"startup_s": round(bare, 4),
            "import_s": round(max(full - bare, 0.0), 4),
            "heavy_import_s": round(max(heavy - bare, 0.0), 4),
            "cold_process_s": round(full, 4),
            "cold_process_heavy_s": round(heavy, 4),
            "module": module}


# ============================================================ the BASELINE arm (s9)
def _baseline_target(p, folds, cfg: P.Config, cache_dir, key):
    """One target through the `s9/final.py` path, with a clock between its own calls.

    This is the reference implementation's work, decomposed -- `s9.final.library_members`,
    `s9.final.windows_all`, its full K x K pairwise matrix, its `synthesise`, its AMBER
    block -- including the `savez_compressed` / `np.load` round trip it really performs
    between `pool`, `synth` and `amber`.  Nothing is reimplemented and nothing is skipped.
    """
    import s9.final as F
    from s7 import audit, debias
    from s7.poolsize import GRID
    from s8.inband import sel_of
    import protein_geometry as geo
    import torsion_lib2 as tl2
    import amber_refine as ar

    clk = P.Clock()
    path = os.path.join(cache_dir, f"{p.pdb}.json")
    with clk("io"):
        done = P._read_json(path)
    if done is not None and done.get("cfg_key") == key:
        done["resumed"] = True
        return done

    fold = folds[p.seq]
    # ---------------- stage 1: retrieve (s9.final.build_pool, first half)
    with clk("windows"):
        peps, frags = F.library_members(p.seq, fold, cfg.n_folds)
        W, PHI, PSI, S, src = F.windows_all(peps + frags, p.n)
    with clk("retrieval"):
        sim = audit.B62[S, audit.encode(p.seq)[None, :]].sum(1)
        order = np.argsort(-sim, kind="stable")
        idx = order[:cfg.k]
        Wk = W[idx]
    # ---------------- stage 2a: the shipped distogram score
    with clk("model_load"):
        F.guard_esm()
        if fold not in F._MODELS:
            import distogram as _dgm
            F._MODELS[fold] = _dgm.train_fold(fold, True, F.N_FOLDS, fragments=True,
                                              verbose=False)
    with clk("distogram"):
        i, j = audit.pair_index(p.n, cfg.min_sep)
        D = audit.pair_dists(Wk, i, j)
        risk = F.distogram_risk(p.seq, fold)
        sc = debias.score_risk(risk, GRID, D)
    # ---------------- stage 2b: the pairwise matrix the filter's medoid consumes
    with clk("pairwise"):
        B = len(Wk)
        Pm = np.zeros((B, B), np.float32)
        for a in range(B):
            Pm[a] = audit.kabsch_rmsd_batch(Wk, Wk[a])
    with clk("reporting"):
        rr = audit.kabsch_rmsd_batch(Wk, p.ca)
        Dnat = audit.pair_dists(np.asarray(p.ca, float)[None], i, j)[0]
    rec = {"pdb": p.pdb, "n": int(p.n), "fold": int(fold), "seq": p.seq,
           "n_windows": int(len(W)), "n_peptides": len(peps), "n_fragments": len(frags),
           "W": Wk.astype(np.float32), "PHI": PHI[idx].astype(np.float32),
           "PSI": PSI[idx].astype(np.float32), "S": S[idx].astype(np.int8),
           "src": src[idx], "sim": sim[idx].astype(np.float32),
           "sc": sc.astype(np.float32), "D": D.astype(np.float32), "P": Pm,
           "rr": rr.astype(np.float32), "Dnat": Dnat.astype(np.float32),
           "nat_ca": np.asarray(p.ca, np.float32)}
    # ---------------- the round trip s9 really performs between its CLI stages
    npz = os.path.join(cache_dir, f"{p.pdb}.npz")
    with clk("serialise"):
        tmp = npz + ".tmp.npz"
        np.savez_compressed(tmp, **{k: np.asarray(v) for k, v in rec.items()})
    with clk("io"):
        os.replace(tmp, npz)
        z = np.load(npz, allow_pickle=True)
        c = {"pdb": str(z["pdb"]), "n": int(z["n"]), "fold": int(z["fold"]),
             "seq": str(z["seq"]), "src": z["src"], "S": z["S"].astype(int),
             "n_windows": int(z["n_windows"])}
        for k in ("W", "PHI", "PSI", "sim", "sc", "D", "P", "rr", "Dnat", "nat_ca"):
            c[k] = z[k].astype(float)
        z.close()
    # ---------------- stage 3: synthesise (filter -> average -> projection)
    v = F.deployable_view(c)
    from s8 import consensus2 as cc
    from s8 import project as pj
    with clk("filter"):
        sub = np.argsort(np.asarray(v["sc"], float), kind="stable")[:cfg.m]
        Ps = np.asarray(v["P"], float)[np.ix_(sub, sub)]
    with clk("average"):
        C = cc._avg_struct(np.asarray(v["W"], float)[sub], Ps)
    with clk("projection"):
        pen = pj.make_penalty(cfg.penalty, v["seq"], int(v["fold"]))
        pth = pj.lam_path(C, pen, (0.0, cfg.lam), maxiter=cfg.maxiter,
                          multi=cfg.multi_start)
        fit, arm = pth[0.0], pth[cfg.lam]
    ca, phi, psi = arm[0], arm[1], arm[2]
    out = {"pdb": p.pdb, "n": c["n"], "fold": c["fold"], "seq": c["seq"],
           "sub": np.asarray(sub, int).tolist(),
           "ca": np.asarray(ca).tolist(), "phi": np.asarray(phi).tolist(),
           "psi": np.asarray(psi).tolist(),
           "fit_ca": np.asarray(fit[0]).tolist(), "avg_ca": np.asarray(C).tolist(),
           "n_windows": c["n_windows"]}
    # ---------------- stage 4: AMBER
    amber_ca = None
    if cfg.amber:
        with clk("amber"):
            BB = geo.build_backbone_batch(np.asarray(phi, float)[None],
                                          np.asarray(psi, float)[None])
            cd = {k: vv[0] for k, vv in BB.items()}
            tab = tl2.library_for(c["seq"], cfg.torsion_window, c["seq"])
            rep = tl2.PerResidueTorsion(c["seq"], tab, chi_bits=False)
            try:
                r = ar.refine_coords(c["seq"], rep, cd, k_restraint=cfg.amber_k,
                                     steps=cfg.amber_steps, components=True)
                amber_ca = np.asarray(r["ca"], float)
                bb = r.get("backbone", {}) or {}
                cm = r.get("components", {}) or {}
                out["amber_e0"] = float(r["energy_initial"])
                out["amber_e1"] = float(r["energy"])
                out["amber_moved"] = float(r["restraint_rmsd"])
                out["amber_strain_after"] = (float(cm.get("bond", 0.0)
                                                   + cm.get("angle", 0.0)) if cm else None)
                out["amber_ca"] = amber_ca.tolist()
                if all(k in bb for k in ("N", "CA", "C")):
                    ph1, ps1 = geo.extract_torsions(bb["N"], bb["CA"], bb["C"])
                    out["amber_phi"] = np.asarray(ph1).tolist()
                    out["amber_psi"] = np.asarray(ps1).tolist()
            except Exception as exc:                                     # noqa: BLE001
                out["amber_err"] = str(exc)[:300]
            ar.clear_cache()
    # ---------------- natives enter ONLY here
    with clk("reporting"):
        nat = c["nat_ca"]
        rrl = c["rr"]
        out["shipped"] = float(sel_of(c["sc"], rrl))
        out["pool_best"] = float(rrl.min())
        out["pool_mean"] = float(rrl.mean())
        out["top_m_best"] = float(rrl[sub].min())
        out["top_m_mean"] = float(rrl[sub].mean())
        out["rmsd_avg"] = float(audit.kabsch_rmsd_batch(np.asarray(C)[None], nat)[0])
        out["rmsd_fit"] = float(audit.kabsch_rmsd_batch(np.asarray(fit[0])[None], nat)[0])
        out["rmsd_arm"] = float(audit.kabsch_rmsd_batch(np.asarray(ca)[None], nat)[0])
        if amber_ca is not None:
            out["rmsd_full"] = float(audit.kabsch_rmsd_batch(amber_ca[None], nat)[0])
            out["d_rmsd"] = out["rmsd_full"] - out["rmsd_arm"]
        else:
            out["rmsd_full"] = None
            out["d_rmsd"] = None
    out["cfg_key"] = key
    out["dev_mode"] = cfg.dev_mode
    with clk("serialise"):
        out["timings"] = {k: round(vv, 4) for k, vv in clk.t.items()}
        out["wall"] = round(clk.total(), 4)
        blob = json.dumps(out)
    with clk("io"):
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            fh.write(blob)
        os.replace(tmp, path)
    out["timings"] = {k: round(vv, 4) for k, vv in clk.t.items()}
    out["wall"] = round(clk.total(), 4)
    out["resumed"] = False
    return out


def run_baseline(manifest_name, cfg: P.Config, verbose=True, limit=None):
    """The `s9/final.py` path, serial, one target at a time -- as it really runs."""
    import s9.final as F
    import peptide_db as db
    tg = P.manifest(manifest_name)
    if limit:
        tg = tg[:limit]
    folds = db.folds(cfg.n_folds)
    tg = sorted(tg, key=lambda p: (folds[p.seq], p.n, p.pdb))
    #: `s9.final.guard_esm` reads `targets()` to decide which missing hot-cache sequence
    #: is a hard error rather than a silent zero vector.  Point it at THIS manifest.
    F.targets = lambda: tg
    F._GUARDED[0] = False
    #: `clock` versions the STAGE ATTRIBUTION, not the science.  Splitting window
    #: extraction out of retrieval, the model load out of scoring and the pairwise matrix
    #: out of filtering changed which bucket a second lands in, and a checkpoint written
    #: under the old attribution cannot be summed with one written under the new.  Bumping
    #: it forces a recompute rather than silently mixing two clocks.
    key = cfg.key({"arm": "baseline_s9", "clock": 2})
    cdir = os.path.join(BASE_CACHE, key)
    os.makedirs(cdir, exist_ok=True)
    rows, errs = [], []
    t0 = time.perf_counter()
    for k, p in enumerate(tg):
        if not P.wait_for_memory(cfg.stop_pct, where=p.pdb):
            errs.append({"pdb": p.pdb,
                         "error": f"memory gate {P.mem_pct()}% for 300s"})
            break
        try:
            r = _baseline_target(p, folds, cfg, cdir, key)
        except Exception:                                                # noqa: BLE001
            import traceback
            errs.append({"pdb": p.pdb, "error": traceback.format_exc()[-1500:]})
            if verbose:
                print(f"  {p.pdb}: FAILED", flush=True)
            continue
        rows.append(r)
        if verbose:
            P._progress(k + 1, len(tg), r)
    return {"manifest": manifest_name, "n_requested": len(tg), "n_done": len(rows),
            "wall": time.perf_counter() - t0, "workers": 1,
            "threads_per_worker": cfg.threads_per_worker, "cfg": cfg.science(),
            "cfg_key": key, "dev_mode": cfg.dev_mode,
            "backends": {"arm": "baseline_s9"}, "cache_dir": cdir, "rows": rows,
            "errors": errs, "stage_totals": P.stage_totals(rows),
            "summary": P.summarise(rows)}


# ============================================================ the run record
#: THIS MACHINE, measured with the AMBER workload itself rather than a proxy
#: (`verify/amber_platform.json`, one process per core on an idle box).  It is an Intel
#: Core Ultra 7 256V: 8 physical cores, NO hyperthreading, and a clean 4/4 P-core /
#: LP-E-core split at 1.67x.
#:
#: WHY EVERY RESULTS FILE CARRIES THIS.  A CPU-second is not a unit of work here: a
#: second on core 6 buys 0.60 of a second on core 1.  So a table showing the optimised
#: arm consuming MORE CPU-seconds than the serial reference is not showing waste, it is
#: showing the same arithmetic done on slower silicon -- if work spreads evenly over all
#: eight cores the predicted inflation is mean(1/rel) = 1.291, and AMBER's measured
#: inflation at 8 workers is 1660.37/1285.85 = 1.291.  Agreement to 0.03%: there is no
#: contention left in that number, only placement.  It also means the parallel ceiling on
#: this box is 6.43x, not 8x, and a ceiling computed as CPU-work/8 is 24% optimistic.
MACHINE = {
    "cpu": "Intel Core Ultra 7 256V",
    "physical_cores": 8, "logical_cores": 8, "hyperthreading": False,
    "relative_throughput_per_core": [0.974, 1.000, 0.935, 0.882,
                                     0.640, 0.693, 0.600, 0.702],
    "core_equivalents": 6.426,
    "fast_cores": 4, "fast_core_equivalents": 3.791,
    #: THE CPU-SECOND INFLATION FLOOR IS SCHEDULE-DEPENDENT, and there are two regimes.
    #: With N cores of unequal speed, the per-structure CPU average depends on how
    #: structures are allotted:
    #:   equal COUNT per core  -> mean(1/rel) = 1.2908   (static equal-count chunks)
    #:   equal WALL per core   -> N/sum(rel) = 1.2449    (ANY dynamic dispatch)
    #: Under dynamic dispatch every core stays busy the same wall time, so fast cores
    #: complete MORE structures and the structure-weighted average shifts toward them.
    #: This is a property of the core map, not of the ordering: a plain work queue gets
    #: it too.  Measured, by two schedules that disagree with each other and each land on
    #: their own regime: contiguous equal-count chunks 1.2913 against 1.2908 (+0.03%),
    #: dynamic longest-first 1.2533 against 1.2449 (+0.67%).  Confirming the core map
    #: twice from two regimes is a stronger result than either agreement alone.
    "inflation_floor_static_equal_count": 1.2908,
    "inflation_floor_dynamic": 1.2449,
    "measured_with": "AMBER refine_coords, production recipe, one process per core",
    "source": "verify/amber_platform.json",
}


def _git_head():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=_ROOT, check=True,
                              capture_output=True, text=True).stdout.strip()[:12]
    except Exception:                                                    # noqa: BLE001
        return None


def _env_report(cfg: P.Config, workers):
    """The threading configuration EVERY timing in this file was taken under.

    `thread_pin` is what the workers were actually started with -- it is set by the parent
    before the pool exists, because OpenBLAS reads its thread count once, at load, and a
    pin applied inside the worker initializer arrives too late.  `thread_env_ambient` is
    what happened to be in the parent's environment; it is recorded for completeness and
    is NOT what governed the run, so nothing should be asserted about it.
    """
    return {"python": sys.version.split()[0], "numpy": np.__version__,
            "platform": sys.platform, "cpu_count": os.cpu_count(),
            "workers": workers, "threads_per_worker": cfg.threads_per_worker,
            "thread_pin": {k: str(cfg.threads_per_worker) for k in P.THREAD_VARS},
            "openmm_threads": 1,
            "openmm_threads_note": (
                "pinned to 1 by amber_refine for DETERMINISM, not for speed: the "
                "hydrogen-frame calibration inherits the thread count and settles in a "
                "different minimum, moving single-point energies by 144 kcal/mol"),
            "thread_env_ambient": {k: os.environ.get(k) for k in P.THREAD_VARS},
            "machine": MACHINE,
            "git": _git_head()}


def _occupancy(rows, workers, wall, cpu_sum):
    """How much of the machine the run actually used, and where the idle went.

    THE AUTOMATIC INVALIDATION SIGNAL.  A run whose workers are idle is either badly
    SCHEDULED or sharing the box, and those are different failures with the same
    symptom in a wall clock.  They separate cleanly once per-worker busy time is known:

        worker-seconds available   = workers x wall
        idle from IMBALANCE        = workers x max_busy - sum(busy)
                                     (what a perfect schedule would have recovered)
        idle from OUTSIDE          = workers x (wall - max_busy)
                                     (process startup, and anything else on the box)

    A run at 62.6% occupancy with a SMALL spread was not badly scheduled; it was sharing
    the machine.  That is exactly the shape of the 445.4 s run this sprint discarded --
    8 workers, 67.4% system CPU, 4 targets whose AMBER correctly declined at 96% RAM --
    and it is flagged here rather than left to be noticed, because the previous version
    of that mistake cost two rounds of suspecting the code.
    """
    per = {}
    for r in rows:
        pid = (r.get("worker") or {}).get("pid")
        if pid is not None:
            per[pid] = per.get(pid, 0.0) + float(r.get("wall") or 0.0)
    if not per or wall <= 0 or workers <= 0:
        return None
    bv = list(per.values())
    avail = workers * wall
    imbalance = max(workers * max(bv) - sum(bv), 0.0)
    outside = max(workers * (wall - max(bv)), 0.0)
    occ = cpu_sum / avail
    spread = (max(bv) - min(bv)) / max(bv) if max(bv) > 0 else 0.0
    return {"effective_cores": round(cpu_sum / wall, 3),
            "occupancy": round(occ, 4),
            "worker_seconds_available": round(avail, 1),
            "worker_seconds_busy": round(sum(bv), 1),
            "idle_from_imbalance_s": round(imbalance, 1),
            "idle_from_outside_s": round(outside, 1),
            "busy_spread_frac": round(spread, 4),
            #: low occupancy WITH a tight spread means the box was shared, not that the
            #: schedule was bad -- the workers were busy, the wall was stretched
            "contention_suspected": bool(occ < 0.75 and spread < 0.15),
            "imbalance_suspected": bool(spread >= 0.15)}


def _worker_footprint(rows):
    """Per-WORKER resident memory, aggregated over the targets each process ran.

    The system-wide RAM trace cannot separate a worker's private pages from the page
    cache every worker shares, and the worker count on this box is decided by the private
    part.  So each target carries its own process's working set and this rolls it up:
    `peak_rss_mb` per pid is what actually has to fit, W times over.
    """
    per = {}
    for r in rows:
        w = r.get("worker") or {}
        pid = w.get("pid")
        if pid is None:
            continue
        d = per.setdefault(pid, {"n_targets": 0, "peak_rss_mb": 0.0, "rss_mb": 0.0,
                                 "window_banks_mb": 0.0})
        d["n_targets"] += 1
        for k in ("peak_rss_mb", "rss_mb", "window_banks_mb"):
            d[k] = max(d[k], float(w.get(k) or 0.0))
    if not per:
        return None
    peaks = [d["peak_rss_mb"] for d in per.values()]
    busy = {pid: 0.0 for pid in per}
    for r in rows:
        pid = (r.get("worker") or {}).get("pid")
        if pid is not None:
            busy[pid] += float(r.get("wall") or 0.0)
    hits = sum((r.get("worker") or {}).get("lib_cache", {}).get("window_hits", 0)
               for r in rows)
    miss = sum((r.get("worker") or {}).get("lib_cache", {}).get("window_misses", 0)
               for r in rows)
    bv = list(busy.values())
    return {"n_processes": len(per), "per_pid": per,
            "busy_s_per_pid": {str(k): round(v, 1) for k, v in busy.items()},
            "busy_s_max": round(max(bv), 1) if bv else None,
            "busy_s_min": round(min(bv), 1) if bv else None,
            "busy_spread_frac": (round((max(bv) - min(bv)) / max(bv), 4)
                                 if bv and max(bv) > 0 else None),
            "peak_rss_mb_max": max(peaks), "peak_rss_mb_mean": float(np.mean(peaks)),
            "peak_rss_gb_total_if_all_resident": sum(peaks) / 1024.0,
            "window_bank_cache": {"hits": hits, "misses": miss,
                                  "hit_rate": (hits / (hits + miss)) if hits + miss else None}}


def _strip(rows):
    """Per-target rows without the coordinate payloads -- the results file stays readable."""
    drop = {"ca", "phi", "psi", "fit_ca", "avg_ca", "amber_ca", "amber_phi", "amber_psi",
            "sub", "seq"}
    return [{k: v for k, v in r.items() if k not in drop} for r in rows]


def bench(arm="optimised", manifest_name="tuning126", workers=1, threads=1,
          fresh=False, limit=None, no_amber=False, out=None, verbose=True,
          backends_env=None, components=False, verify_arms=False,
          amber_concurrency=0):
    if components and arm == "baseline":
        raise SystemExit("the baseline arm IS `s9/final.py`, which contains no VQE and "
                         "no Legacy stage; --components is only meaningful on the "
                         "optimised arm, and the two end-to-end numbers must be "
                         "reported separately.")
    cfg = replace(P.PROD, threads_per_worker=threads, amber=not no_amber,
                  quantum=components, legacy=components,
                  report_single_start_fit=verify_arms,
                  amber_concurrency=amber_concurrency,
                  dev_mode=("no-amber" if no_amber else P.PROD.dev_mode))
    if manifest_name.startswith("smoke"):
        cfg = replace(cfg, dev_mode=(cfg.dev_mode + "+" if cfg.dev_mode else "")
                      + manifest_name)
    if fresh:
        _clear_cache(arm, cfg, backends_env)

    startup = measure_startup("s9.final" if arm == "baseline" else "core.pipeline")
    t_wall0 = time.perf_counter()
    with Resources() as res:
        if arm == "baseline":
            out_run = run_baseline(manifest_name, cfg, verbose=verbose, limit=limit)
        else:
            out_run = P.run(manifest_name, cfg, workers=workers, verbose=verbose,
                            backends_env=backends_env, limit=limit)
    wall = time.perf_counter() - t_wall0

    rows = out_run.pop("rows")
    st = out_run["stage_totals"]
    rec = {
        "arm": arm,
        "manifest": manifest_name,
        "n_requested": out_run["n_requested"],
        "n_done": out_run["n_done"],
        "n_fresh": st["_n_fresh"],
        "n_resumed": st["_n_resumed"],
        "cold": bool(st["_n_resumed"] == 0),
        "dev_mode": cfg.dev_mode,
        "headline_eligible": bool(
            cfg.dev_mode == "" and st["_n_resumed"] == 0
            and out_run["n_done"] == out_run["n_requested"]
            #: and the box was not shared.  A wall clock measured while another process
            #: had the machine is not this pipeline's wall clock.
            and not (_occupancy(rows, out_run["workers"], wall, st["_sum"]) or {}
                     ).get("contention_suspected", False)
            #: and every target must have produced every arm.  A target whose AMBER
            #: declined at the memory ceiling still counts in `n_done` -- the failure is
            #: recorded per target, not as a run error -- so a run can be "complete" and
            #: still be missing an arm on a fifth of the manifest.  That moves an arm's
            #: MEAN without moving any number, which is the subtlest way this harness
            #: could mislead, so it disqualifies the run.
            and not any(t.get("amber_err") for t in rows)),
        "n_amber_declined": sum(1 for t in rows if t.get("amber_err")),
        #: which pipeline this number is FOR.  "like-for-like" is the s9/final.py path and
        #: is the only arm a speedup against s9/final.py may be quoted from; "four-
        #: component" adds the VQE selector and the Legacy refiner, which are work the
        #: reference does not do, so its wall clock is a different experiment's.
        "system": ("four-component" if components
                   else ("like-for-like-s9+verification" if verify_arms
                         else "like-for-like-s9")),
        "verification_arms": bool(verify_arms),
        "components": {"amber": cfg.amber, "quantum": cfg.quantum,
                       "legacy": cfg.legacy},
        "wall_s": round(wall, 3),
        "wall_run_s": round(out_run["wall"], 3),
        "workers": out_run["workers"],
        "threads_per_worker": cfg.threads_per_worker,
        "process_costs": startup,
        #: per-stage seconds over the WHOLE manifest, each measured when its target was
        #: computed.  Valid whether or not this invocation had to resume.
        "stages_s": {k: round(st[k], 4) for k in P.STAGES},
        "stage_cpu_sum_s": round(st["_sum"], 3),
        #: the same, restricted to what this invocation actually computed
        "stages_fresh_s": {k: round(st[f"_fresh_{k}"], 4) for k in P.STAGES},
        "stage_cpu_sum_fresh_s": round(st["_fresh_sum"], 3),
        "per_target_wall_s": (round(wall / out_run["n_done"], 4)
                              if out_run["n_done"] else None),
        "throughput_targets_per_min": (round(60.0 * out_run["n_done"] / wall, 3)
                                       if wall > 0 else None),
        "resources": res.report(),
        "worker_footprint": _worker_footprint(rows),
        "occupancy": _occupancy(rows, out_run["workers"], wall, st["_sum"]),
        "config": out_run["cfg"],
        "cfg_key": out_run["cfg_key"],
        "backends": out_run["backends"],
        "cache_dir": out_run["cache_dir"],
        "science": out_run["summary"],
        "errors": out_run["errors"],
        "env": _env_report(cfg, out_run["workers"]),
        "per_target": _strip(rows),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    os.makedirs(RESULTS, exist_ok=True)
    name = out or f"{arm}_{manifest_name}" + (f"_w{workers}" if arm != "baseline" else "")
    path = os.path.join(RESULTS, name + ".json")
    P._atomic_write_json(path, rec)
    if verbose:
        print_run(rec)
        print(f"\nwritten: {path}", flush=True)
    return rec, path


def _clear_cache(arm, cfg, backends_env):
    import shutil
    if arm == "baseline":
        d = os.path.join(BASE_CACHE, cfg.key({"arm": "baseline_s9", "clock": 2}))
    else:
        d = os.path.join(P.CACHE_ROOT,
                         cfg.key(core.backend_report() if backends_env is None else None))
    if os.path.isdir(d):
        shutil.rmtree(d, ignore_errors=True)


# ============================================================ printing
#: (key, label).  The first block is the like-for-like `s9/final.py` pipeline; the second
#: is the two mandated components that file does not contain.
ARM_NAMES = (
    ("shipped", "shipped argmin"),
    ("pool_best", "pool best (ceiling)"),
    ("top_m_best", "top-75 best (ceiling)"),
    ("rmsd_avg", "raw average (illegal)"),
    ("rmsd_fit", "projection, no prior"),
    ("rmsd_fit_single", "  single-start fit"),
    ("rmsd_arm", "SYNTHESIS (stage 3)"),
    ("rmsd_full", "S9 SYSTEM (+AMBER)"),
    ("rmsd_vqe_sel_uniform", "  uniform ensemble sel"),
    ("rmsd_vqe_sel", "  CVaR-VQE selection"),
    ("rmsd_legacy_centre", "  neighbourhood centre"),
    ("rmsd_legacy_sel", "  Legacy late refiner"),
    ("rmsd_u_synth", "  uniform 2^n synthesis"),
    ("rmsd_q_synth", "  quantum synthesis"),
    ("rmsd_system", "FOUR-COMPONENT SYSTEM"),
)


def print_run(r):
    print("\n" + "=" * 78)
    print(f"{r['arm'].upper()}  manifest={r['manifest']}  n={r['n_done']}/"
          f"{r['n_requested']}  workers={r['workers']} x "
          f"{r['threads_per_worker']} thread(s)")
    print(f"system: {r.get('system', 'like-for-like-s9')}   components "
          f"{r.get('components', {})}")
    print("=" * 78)
    if not r["headline_eligible"]:
        why = []
        if r["dev_mode"]:
            why.append(f"dev_mode={r['dev_mode']}")
        if r["n_resumed"]:
            why.append(f"{r['n_resumed']} targets resumed from cache")
        if r["n_done"] != r["n_requested"]:
            why.append(f"{r['n_requested'] - r['n_done']} targets missing")
        if r.get("n_amber_declined"):
            why.append(f"AMBER declined on {r['n_amber_declined']} targets "
                       f"(memory ceiling) -- an arm is missing, its mean is over a subset")
        oc = r.get("occupancy") or {}
        if oc.get("contention_suspected"):
            why.append(f"occupancy {oc['occupancy']*100:.0f}% with a "
                       f"{oc['busy_spread_frac']*100:.0f}% worker spread: the box was shared")
        print(f"NOT A HEADLINE NUMBER: {'; '.join(why)}")
    tot = sum(r["stages_s"].values()) or 1.0
    print(f"\n{'stage':14}{'CPU-s':>10}{'share':>8}   (summed over {r['n_done']} targets, "
          f"{r['n_fresh']} computed now)")
    for k, v in r["stages_s"].items():
        print(f"{k:14}{v:10.2f}{100*v/tot:7.1f}%")
    print(f"{'-'*14}{'-'*10}")
    print(f"{'stage sum':14}{tot:10.2f}")
    pc = r["process_costs"]
    nproc = max(r["workers"], 1)
    print(f"{'startup':14}{pc['startup_s']:10.2f}   (per process, x{nproc})")
    print(f"{'imports':14}{pc['import_s']:10.2f}   ({pc['module']}, per process)")
    print(f"{'imports(heavy)':14}{pc['heavy_import_s']:10.2f}   "
          f"(+torch/openmm/db, per process)")
    print(f"{'process tax':14}{nproc * pc['cold_process_heavy_s']:10.2f}   "
          f"(x{nproc} workers, cold)")
    print(f"\nwall {r['wall_s']:.1f}s   per target {r['per_target_wall_s']}s   "
          f"throughput {r['throughput_targets_per_min']}/min")
    res = r["resources"]
    if res.get("cpu_pct"):
        print(f"CPU {res['cpu_pct']['mean']:.1f}% mean / {res['cpu_pct']['peak']:.1f}% peak"
              f"   RAM {res['ram_pct']['mean']:.1f}% mean / {res['ram_pct']['peak']:.1f}% "
              f"peak ({res['ram_used_gb']['peak']:.2f} GB of "
              f"{res['total_ram_gb']:.1f} GB)   n={res['n_samples']} samples")
    s = r["science"]
    if s:
        print("\narm                       mean   median       sd      min      max   <2.0")
        for k, name in ARM_NAMES:
            d = s.get(k)
            if not d:
                continue
            print(f"{name:24}{d['mean']:8.4f}{d['median']:9.4f}{d['sd']:9.4f}"
                  f"{d['min']:9.3f}{d['max']:9.3f}{d['frac_under_2.0']:7.3f}")
        e = s.get("vqe_entropy_bits")
        if e:
            print(f"\nCVaR state entropy: mean {e['mean']:.2f} bits, min "
                  f"{e['min']:.2f}; {e['n_collapsed']} target(s) collapsed below 0.5 bits "
                  f"(a collapsed state IS the shipped argmin and contributes nothing)")
        dup = s.get("duplicate_candidates")
        if dup:
            print(f"duplicate candidates in the filtered prefix: "
                  f"{dup['frac_duplicate']*100:.1f}% "
                  f"({dup['n_distinct_mean']:.0f} distinct of "
                  f"{dup['n_top_mean']:.0f}); the pairwise stage computes one row per "
                  f"distinct structure")
        ab = s.get("ablations") or {}
        if ab:
            print("\n--- component ablations (negative = the component helps) ---")
            print(f"{'ablation':32}{'with':>8}{'without':>9}{'d':>9}"
                  f"{'95% CI':>20}{'W/L':>9}")
            for name, v in ab.items():
                print(f"{name:32}{v['mean_with']:8.4f}{v['mean_without']:9.4f}"
                      f"{v['mean_diff']:+9.4f}"
                      f"  [{v['ci95'][0]:+.4f},{v['ci95'][1]:+.4f}]"
                      f"{v['n_better']:5d}/{v['n_worse']:<4d}")
                if v["contribution_is_zero_within_ci"]:
                    print(f"{'':32}  contribution is ZERO within its CI -- reported as a "
                          f"result, not hidden")
    oc = r.get("occupancy")
    if oc:
        print(f"occupancy {oc['occupancy']*100:.1f}% ({oc['effective_cores']:.2f} of "
              f"{r['workers']} workers busy)   idle {oc['idle_from_imbalance_s']:.0f} s "
              f"imbalance + {oc['idle_from_outside_s']:.0f} s outside   "
              f"worker spread {oc['busy_spread_frac']*100:.1f}%")
        if oc["contention_suspected"]:
            print("  ** CONTENTION SUSPECTED: occupancy is low and the workers finished "
                  "together.\n     The workers were busy and the WALL was stretched -- "
                  "another process had the box.\n     This run is not headline-eligible "
                  "and its wall clock is not this pipeline's.")
        elif oc["imbalance_suspected"]:
            print("  ** IMBALANCE: the workers did not finish together.  Real, and the "
                  "schedule's to fix.")
    wf = r.get("worker_footprint")
    if wf:
        print(f"\nper-worker footprint: {wf['n_processes']} process(es), peak RSS "
              f"{wf['peak_rss_mb_mean']:.0f} MB mean / {wf['peak_rss_mb_max']:.0f} MB max"
              f"  ({wf['peak_rss_gb_total_if_all_resident']:.2f} GB if all resident)")
        c = wf["window_bank_cache"]
        if c["hit_rate"] is not None:
            print(f"window-bank cache: {c['hits']} hits / {c['misses']} misses "
                  f"({c['hit_rate']*100:.0f}%)")
    if r["errors"]:
        print(f"\n{len(r['errors'])} error(s): "
              f"{[e['pdb'] for e in r['errors']][:8]}")


def compare(a_path, b_path, verbose=True):
    """Before/after/speedup, per stage, from two real runs on disk."""
    a = P._read_json(a_path)
    b = P._read_json(b_path)
    if a is None or b is None:
        raise SystemExit(f"cannot read {a_path if a is None else b_path}")
    if a["manifest"] != b["manifest"]:
        print(f"WARNING: different manifests ({a['manifest']} vs {b['manifest']}); "
              f"a per-stage speedup across different experiments is not a speedup.")
    if a["n_done"] != b["n_done"]:
        print(f"WARNING: {a['n_done']} vs {b['n_done']} targets completed; the per-stage "
              f"sums are over different experiments.")
    if not (a["cold"] and b["cold"]):
        print("NOTE: at least one arm resumed targets from its checkpoints.  The "
              "per-stage table is still valid -- each target's stage clock was taken "
              "when that target was computed -- but the WALL CLOCK row is a property "
              "of one invocation and is NOT a cold end-to-end measurement.")
    rows = []
    for k in P.STAGES:
        x, y = a["stages_s"].get(k, 0.0), b["stages_s"].get(k, 0.0)
        rows.append((k, x, y, (x / y) if y > 1e-9 else (float("inf") if x > 1e-9 else 1.0)))
    out = {"baseline": a_path, "optimised": b_path, "manifest": a["manifest"],
           "n": {"baseline": a["n_done"], "optimised": b["n_done"]},
           "cold": {"baseline": a["cold"], "optimised": b["cold"]},
           "headline_eligible": {"baseline": a["headline_eligible"],
                                 "optimised": b["headline_eligible"]},
           "workers": {"baseline": a["workers"], "optimised": b["workers"]},
           "stage_speedup": {k: {"baseline_s": round(x, 3), "optimised_s": round(y, 3),
                                 "speedup": (round(s, 3) if s != float("inf") else None)}
                             for k, x, y, s in rows},
           "wall_s": {"baseline": a["wall_s"], "optimised": b["wall_s"]},
           "end_to_end_speedup": round(a["wall_s"] / b["wall_s"], 3) if b["wall_s"] else None,
           "cpu_stage_sum_speedup": (round(a["stage_cpu_sum_s"] / b["stage_cpu_sum_s"], 3)
                                     if b["stage_cpu_sum_s"] else None),
           #: CPU-seconds are not comparable across worker counts on this machine (see
           #: MACHINE).  These say how much of the box's REAL capacity each arm used.
           "core_equivalents_available": MACHINE["core_equivalents"],
           "effective_cores": {"baseline": round(a["stage_cpu_sum_s"] / a["wall_s"], 3),
                               "optimised": round(b["stage_cpu_sum_s"] / b["wall_s"], 3)},
           #: `baseline_cpu / core_equivalents`: what perfectly spreading the REFERENCE
           #: workload over this machine would cost.  The optimised arm may sit under it,
           #: because the exact structural wins removed CPU-work before parallelising.
           "machine_ceiling_s": round(a["stage_cpu_sum_s"]
                                      / MACHINE["core_equivalents"], 2),
           "machine_ceiling_speedup": round(a["wall_s"] * MACHINE["core_equivalents"]
                                            / a["stage_cpu_sum_s"], 3),
           "cpu_work_removed_s": round(a["stage_cpu_sum_s"] - b["stage_cpu_sum_s"], 2),
           "science": {"baseline": _sci(a), "optimised": _sci(b)},
           "science_delta": _sci_delta(a, b),
           "resources": {"baseline": a["resources"], "optimised": b["resources"]}}
    if verbose:
        print("\n" + "=" * 78)
        print(f"PER-STAGE  baseline -> optimised   manifest={out['manifest']}  "
              f"n={a['n_done']} vs {b['n_done']}")
        print("=" * 78)
        print(f"{'stage':14}{'baseline s':>12}{'optimised s':>13}{'speedup':>10}")
        for k, x, y, s in rows:
            sv = "  --" if s == float("inf") else f"{s:9.2f}x"
            print(f"{k:14}{x:12.2f}{y:13.2f}{sv:>10}")
        print(f"{'-'*49}")
        print(f"{'stage sum':14}{a['stage_cpu_sum_s']:12.2f}"
              f"{b['stage_cpu_sum_s']:13.2f}{out['cpu_stage_sum_speedup']:9.2f}x")
        print(f"{'WALL CLOCK':14}{a['wall_s']:12.2f}{b['wall_s']:13.2f}"
              f"{out['end_to_end_speedup']:9.2f}x   "
              f"({a['workers']}w -> {b['workers']}w)")
        ce = MACHINE["core_equivalents"]
        eff = out["effective_cores"]
        print(f"\nA CPU-SECOND IS NOT A UNIT OF WORK ON THIS BOX.  {MACHINE['cpu']}, "
              f"{MACHINE['physical_cores']} cores, no hyperthreading,\n"
              f"a 4/4 P-core / LP-E split at 1.67x -- so the machine is {ce} "
              f"fastest-core-equivalents, not {MACHINE['physical_cores']}.\n"
              f"Where the optimised arm shows MORE CPU-seconds than the serial "
              f"reference, that is the same\narithmetic on slower silicon, not waste.")
        print(f"effective cores {eff['baseline']:.2f} -> {eff['optimised']:.2f}")
        #: THE CEILING, and why the achieved number can sit at or slightly past it.
        #: The baseline is one worker, so Windows parks it on a fast core and its
        #: CPU-seconds ARE fast-core-seconds.  Spreading exactly that work over the whole
        #: machine takes `baseline_cpu / core_equivalents`.  The optimised arm beats that
        #: only because it does LESS work: the exact structural wins removed real
        #: CPU-seconds before any of it was parallelised.  Both terms are printed so the
        #: two effects are never conflated into a single efficiency figure.
        ceiling = a["stage_cpu_sum_s"] / ce
        saved = a["stage_cpu_sum_s"] - b["stage_cpu_sum_s"]
        print(f"machine ceiling on the REFERENCE workload: {a['stage_cpu_sum_s']:.0f} "
              f"fast-core-s / {ce} = {ceiling:.1f} s ({a['wall_s'] / ceiling:.2f}x)")
        print(f"achieved {b['wall_s']:.1f} s ({out['end_to_end_speedup']:.2f}x) -- at "
              f"that ceiling, and past it only because the exact structural wins removed "
              f"{saved:.0f} CPU-s\nof work before any of it was parallelised.  "
              f"Scheduling has essentially nothing left; what remains is CPU-work.")
        print("\nscience (must not move):")
        print(f"{'arm':24}{'baseline':>11}{'optimised':>11}{'delta':>11}")
        for k in ("shipped", "pool_best", "rmsd_fit", "rmsd_arm", "rmsd_full"):
            x = out["science"]["baseline"].get(k)
            y = out["science"]["optimised"].get(k)
            if x is None or y is None:
                continue
            print(f"{k:24}{x:11.4f}{y:11.4f}{y - x:+11.6f}")
        for k, v in out["resources"].items():
            c = v.get("cpu_pct") or {}
            m = v.get("ram_pct") or {}
            print(f"{k:10} CPU {c.get('mean', 0):5.1f}% mean {c.get('peak', 0):5.1f}% peak"
                  f"   RAM {m.get('mean', 0):5.1f}% mean {m.get('peak', 0):5.1f}% peak")
    return out


def _sci(r):
    s = r.get("science") or {}
    return {k: v["mean"] for k, v in s.items() if isinstance(v, dict) and "mean" in v}


def _sci_delta(a, b):
    x, y = _sci(a), _sci(b)
    return {k: (y[k] - x[k]) for k in x if k in y}


# ============================================================ worker sweep
def sweep_workers(counts, manifest_name="smoke24", threads=1, no_amber=False,
                  verbose=True, amber_concurrency=0):
    """Find the worker count experimentally.  Each point is a real, cold run.

    The trap this exists to catch: with W workers each spawning T BLAS/OpenMP threads the
    box runs W*T runnable threads on `cpu_count` cores and throughput FALLS while CPU%
    LOOKS better.  So every point reports throughput, CPU% and RAM peak together, and the
    winner is chosen on throughput with RAM peak under the gate -- never on CPU%.
    """
    pts = []
    for w in counts:
        rec, path = bench("optimised", manifest_name, workers=w, threads=threads,
                          fresh=True, no_amber=no_amber, verbose=False,
                          amber_concurrency=amber_concurrency,
                          out=f"sweep_{manifest_name}_w{w}t{threads}"
                              + (f"a{amber_concurrency}" if amber_concurrency else ""))
        pts.append({"workers": w, "threads": threads, "wall_s": rec["wall_s"],
                    "throughput": rec["throughput_targets_per_min"],
                    "cpu_pct_mean": rec["resources"]["cpu_pct"]["mean"],
                    "ram_pct_peak": rec["resources"]["ram_pct"]["peak"],
                    "ram_used_gb_peak": rec["resources"]["ram_used_gb"]["peak"],
                    "stage_cpu_sum_s": rec["stage_cpu_sum_s"],
                    "science": _sci(rec), "path": path})
        if verbose:
            p = pts[-1]
            print(f"workers={w:2d} threads={threads}  wall {p['wall_s']:8.1f}s  "
                  f"{p['throughput']:6.2f} targets/min  CPU {p['cpu_pct_mean']:5.1f}%  "
                  f"RAM peak {p['ram_pct_peak']:5.1f}% ({p['ram_used_gb_peak']:.2f} GB)",
                  flush=True)
    out = sweep_verdict(pts, manifest_name, threads)
    P._atomic_write_json(os.path.join(RESULTS, f"worker_sweep_{manifest_name}.json"), out)
    if verbose:
        _print_sweep(out)
    return out


def sweep_verdict(pts, manifest_name, threads, gate=92):
    """Choose the operating point, and verify the NUMBERS at every point, per target.

    Two separate questions, which a mean over an arm cannot tell apart:

      1. DID ANY NUMBER MOVE?  Compared per target against the 1-worker run, on the
         targets where both produced a value.  Parallelism is a scheduling decision and
         must move nothing; anything but 0.0 here is a hard failure, and it is the check
         that catches the real hazard -- OpenMM's minimisation result depends on its
         thread environment, so a thread count is not only a speed knob.
      2. DID EVERY TARGET FINISH?  A configuration that runs the box out of memory makes
         the AMBER stage decline to start a context, which leaves an arm null and moves
         its MEAN without moving a single number.  That is a coverage failure, not a
         numerical one, and reporting it as "the science changed" would be wrong.

    The chosen point is the fastest one with FULL coverage and a RAM peak under the gate.
    Never the highest CPU%: an oversubscribed box shows more CPU and less work.
    """
    base = _read_point(pts[0]["path"])
    #: only the arms this configuration actually produces.  With `--components` off,
    #: `rmsd_vqe_sel` and friends are null on every target, and counting them as missing
    #: would report every point as incomplete.
    live_arms = [a for a in P.ARMS
                 if any(t.get(a) is not None for t in base.values())]
    rows = []
    for p in pts:
        r = _read_point(p["path"])
        moved, cov = {}, {}
        for arm in live_arms:
            common = [(base[t].get(arm), r[t].get(arm)) for t in r
                      if t in base and base[t].get(arm) is not None
                      and r[t].get(arm) is not None]
            if common:
                moved[arm] = max(abs(x - y) for x, y in common)
            cov[arm] = sum(1 for t in r if r[t].get(arm) is not None)
        row = dict(p)
        row["worst_abs_diff_vs_1_worker"] = (max(moved.values()) if moved else 0.0)
        row["moved_per_arm"] = moved
        row["n_targets_complete"] = min(cov.values()) if cov else 0
        row["n_targets"] = len(r)
        row["full_coverage"] = bool(cov and min(cov.values()) == len(r))
        row["under_memory_gate"] = bool(p["ram_pct_peak"] < gate)
        rows.append(row)
    ok = [x for x in rows if x["full_coverage"] and x["under_memory_gate"]]
    best = max(ok or rows, key=lambda x: x["throughput"] or 0.0)
    return {"manifest": manifest_name, "threads_per_worker": threads,
            "cpu_count": os.cpu_count(), "memory_gate_pct": gate, "points": rows,
            "best_workers": best["workers"], "best_throughput": best["throughput"],
            "best_chosen_from": "full coverage and RAM peak under the gate",
            "numbers_identical_at_every_worker_count": bool(
                all(x["worst_abs_diff_vs_1_worker"] == 0.0 for x in rows)),
            "worst_abs_diff_over_all_points": max(
                x["worst_abs_diff_vs_1_worker"] for x in rows),
            "points_with_incomplete_coverage": [x["workers"] for x in rows
                                                if not x["full_coverage"]],
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")}


def _read_point(path):
    r = P._read_json(path) or {}
    return {t["pdb"]: t for t in r.get("per_target", [])}


def _print_sweep(out):
    print(f"\n{'workers':>8}{'wall s':>10}{'targets/min':>13}{'CPU%':>8}{'RAM peak%':>11}"
          f"{'complete':>10}{'max |d| vs w=1':>16}")
    for p in out["points"]:
        print(f"{p['workers']:8d}{p['wall_s']:10.1f}{p['throughput']:13.2f}"
              f"{p['cpu_pct_mean']:8.1f}{p['ram_pct_peak']:11.1f}"
              f"{p['n_targets_complete']:6d}/{p['n_targets']:<3d}"
              f"{p['worst_abs_diff_vs_1_worker']:16.2e}")
    print(f"\nCHOSEN: {out['best_workers']} workers x {out['threads_per_worker']} "
          f"thread(s) -- {out['best_chosen_from']}")
    print(f"numbers identical at every worker count: "
          f"{out['numbers_identical_at_every_worker_count']} "
          f"(worst |diff| {out['worst_abs_diff_over_all_points']:.2e})")
    inc = out["points_with_incomplete_coverage"]
    if inc:
        print(f"incomplete coverage at {inc} worker(s): the memory gate stopped the "
              f"AMBER stage from starting a context.  That moves an arm's MEAN without "
              f"moving any number, and is a coverage failure, not a numerical one.")


# ============================================================ cli
def main(argv=None):
    ap = argparse.ArgumentParser(prog="core.bench",
                                 description=__doc__.split("\n")[1])
    ap.add_argument("--arm", choices=("baseline", "optimised"), default="optimised")
    ap.add_argument("--manifest", default="tuning126", choices=P.MANIFESTS)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--fresh", action="store_true",
                    help="delete this arm's checkpoints first: a cold measurement")
    ap.add_argument("--no-amber", action="store_true",
                    help="DEV MODE: skip the validity stage; excluded from headlines")
    ap.add_argument("--amber-concurrency", type=int, default=0,
                    help="cap how many AMBER minimisations run at once across all "
                         "workers (0 = no cap).  Admission control only: it changes WHEN "
                         "work runs, never what is computed.")
    ap.add_argument("--verify-arms", action="store_true",
                    help="also emit the single-start projection, so BOTH published "
                         "incumbents (3.2005 and 3.204) reproduce from one run.  It is "
                         "work s9/final.py does not do, so a run with it on is a "
                         "VERIFICATION pass and not a headline timing run.")
    ap.add_argument("--components", action="store_true",
                    help="the full four-component system (adds the CVaR-VQE selector and "
                         "the Legacy refiner).  This is MORE work than s9/final.py does, "
                         "so its wall clock is a different experiment's and is labelled "
                         "as such in the results file.")
    ap.add_argument("--out", default=None)
    ap.add_argument("--backends", default=None)
    ap.add_argument("--compare", nargs=2, metavar=("BASELINE", "OPTIMISED"))
    ap.add_argument("--sweep-workers", default=None,
                    help="comma-separated worker counts to measure, e.g. 1,2,4,6,8")
    ap.add_argument("--i-am-spending-the-benchmark", action="store_true")
    a = ap.parse_args(argv)

    if a.compare:
        compare(a.compare[0], a.compare[1])
        return 0
    if a.manifest == "benchmark60" and not a.i_am_spending_the_benchmark:
        print("REFUSED: the 60-target benchmark's single pre-registered pass is spent. "
              "Optimisation work is measured on tuning126.")
        return 2
    if a.sweep_workers:
        sweep_workers([int(x) for x in a.sweep_workers.split(",")], a.manifest,
                      threads=a.threads, no_amber=a.no_amber,
                      amber_concurrency=a.amber_concurrency)
        return 0
    bench(a.arm, a.manifest, workers=a.workers, threads=a.threads, fresh=a.fresh,
          limit=a.limit, no_amber=a.no_amber, out=a.out, backends_env=a.backends,
          components=a.components, verify_arms=a.verify_arms,
          amber_concurrency=a.amber_concurrency)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
