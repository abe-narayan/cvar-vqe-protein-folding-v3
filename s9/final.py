"""Sprint 9 FINAL EVALUATION -- the one clean pass on the 60-target held-out benchmark.

THE PRE-REGISTRATION
====================
Everything in this file was fixed before a single benchmark number was read.  Nothing here
may be tuned on these 60 targets; every hyperparameter is inherited from the 126-target
tuning instrument and the pre-registration restates it:

    K = 500                the BLOSUM62 retrieval depth  (S7, S8-5)
    filter = distogram score's top 75                    (S8-11)
    synthesis = coordinate average -> projection onto the manifold of ideal-geometry
                chains, MULTI-START at every rung, with the hinged Ramachandran prior
                `ramah` at lambda = 0.3                  (S8-11, S9-2)
    validity  = AMBER ff14SB/GBn2 restrained relaxation, k = 10 kcal/mol/A^2,
                steps = 0 (converged)                    (S8-12, S9-5)

THE FOUR STAGES, and what each is for
-------------------------------------
1. RETRIEVE   Out-of-fold peptides + `distogram._fold_fragments(fold, 5)` are windowed at
              the target length; the K=500 highest BLOSUM62 sums are the pool.  Identical
              construction to `s7/audit.py:build_target` (the reference implementation),
              including the *stable* argsort, because the similarity has large tie sets
              and a different tie-break gives a different pool.
2. FILTER     The shipped distogram Bayes-risk score's top 75.  S8-11 measured that this
              filter's job is collapsing the pool onto one mode, not keeping good members.
3. SYNTHESISE Coordinate-average the 75 after superposing them on their medoid, then
              project that average back onto the manifold of ideal-geometry chains.  The
              average has a mean CA-CA bond of ~2.96 A and is not a peptide; the
              projection is the correct inverse (rescaling is 1.04 A worse).  Multi-start
              is load-bearing: the projection is DEGENERATE -- a CA trace admits two
              ideal-geometry torsion solutions at near-equal distance, one plausible and
              one not -- and warm-starting cannot cross between them (S9-2).
4. RELAX      One restrained ff14SB/GBn2 minimisation.  This is a VALIDITY stage and it
              COSTS accuracy (+0.011 to +0.026 A); it earns its place by removing 10^4 to
              10^6 kcal/mol of builder strain so that nothing downstream inherits it.

The PRIMARY reported number is the full four-stage system INCLUDING AMBER.  Stage 3 alone
is reported alongside as the declared ablation.  Both go in the record.

WHAT IS A BASELINE AND WHAT IS A DIAGNOSTIC
-------------------------------------------
Deployable, reported as baselines:  the shipped distogram argmin (the historical
pipeline), and the unfiltered pool's own best member (the achievability CEILING -- not
deployable, but the number that prices the selection gap).  The perfect-distance-oracle
ranks the pool by L1 agreement with the NATIVE distance matrix; it is a diagnostic and is
labelled as one everywhere it appears.

NATIVES
-------
`nat_ca`, `rr` and `Dnat` are REPORTING LABELS.  Nothing on the deployable path reads
them: `deployable_view()` is the only object the synthesis and relaxation see, and
`stage_leak` NaN-poisons every native quantity and asserts the emitted structures are
bit-identical.

STAGES (each resumable, each writing atomically)
    python -m s9.final pool     # K=500 pools + distogram scores + pairwise RMSD
    python -m s9.final synth    # filter 75 -> average -> projection (ramah@0.3, multi)
    python -m s9.final amber    # the restrained relaxation, one context at a time
    python -m s9.final audit    # leakage: identities, intersections, fold discipline
    python -m s9.final leak     # NaN-poisoning identity test
    python -m s9.final report   # every table in the mandate
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np                                                      # noqa: E402

os.environ.setdefault("NT", "2")
try:
    import torch                                                        # noqa: E402
    torch.set_num_threads(2)
except Exception:                                                       # noqa: BLE001
    pass

import distogram as dgm                                                 # noqa: E402
import peptide_db as db                                                 # noqa: E402
import protein_geometry as geo                                          # noqa: E402
from s7 import audit, debias                                            # noqa: E402
from s7.poolsize import GRID                                            # noqa: E402
from s8 import consensus2 as cc                                         # noqa: E402
from s8 import project as pj                                            # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "final_cache")

# ---------------------------------------------------------------- PRE-REGISTERED CONSTANTS
K = 500               #: retrieval depth
M = 75                #: distogram-score filter size
PEN = "ramah"         #: the hinged Ramachandran prior
LAM = 0.3             #: its weight
AMBER_K = 10.0        #: restraint force constant, kcal/mol/A^2
AMBER_STEPS = 0       #: 0 = minimise to convergence (OpenMM convention)
N_FOLDS = 5

#: the two numbers this pass is asked to replicate
TUNING_MEAN = 3.204   #: 126-target tuning instrument
DEV_MEAN = 3.221      #: 24-target dev set
#: and the EFFECT behind them -- S8-11's synthesis against the shipped distogram argmin.
#: The absolute mean is a property of how hard a target set is; the paired gain is the
#: property of the ARCHITECTURE, so this is what "did it replicate" actually asks about.
TUNING_GAIN = -0.250              #: 126 targets, [-0.383, -0.117], 80W/46L
TUNING_GAIN_CI = (-0.383, -0.117)
DEV_GAIN = -0.255                 #: 24 targets, [-0.501, -0.008], 16W/8L
TUNING_BASELINE = 3.454           #: the shipped distogram argmin on the tuning instrument
TUNING_POOL_BEST = 1.711          #: its achievability ceiling

#: resource discipline -- see the mandate.  ctypes GlobalMemoryStatusEx, not a subprocess.
START_PCT = int(os.environ.get("START_PCT", "78"))
STOP_PCT = int(os.environ.get("STOP_PCT", "90"))


# ============================================================ resource discipline
def mem_pct():
    """Percent of physical RAM in use, from a ctypes GlobalMemoryStatusEx read."""
    import ctypes

    class _MS(ctypes.Structure):
        _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    st = _MS()
    st.dwLength = ctypes.sizeof(_MS)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st))
    return int(st.dwMemoryLoad)


def drop_caches():
    import gc
    try:
        import amber_refine as ar
        ar.clear_cache()
    except Exception:                                                   # noqa: BLE001
        pass
    gc.collect()


def start_gate():
    m = mem_pct()
    if m >= START_PCT:
        print(f"MEMORY {m}% >= {START_PCT}%: refusing to start. Every stage is resumable "
              f"per target; rerun when the box is quieter.", flush=True)
        return False
    print(f"memory {m}% at start", flush=True)
    return True


def gate(where=""):
    m = mem_pct()
    if m < STOP_PCT:
        return True
    drop_caches()
    m2 = mem_pct()
    print(f"MEMORY {m}% -> {m2}% after dropping caches; "
          f"{'continuing' if m2 < STOP_PCT else 'STOPPING'} {where}", flush=True)
    return m2 < STOP_PCT


# ============================================================ small io
def _p(name):
    return os.path.join(HERE, name)


def _write(name, obj):
    tmp = _p(name) + ".tmp"
    for attempt in range(5):
        try:
            with open(tmp, "w") as fh:
                json.dump(obj, fh)
            os.replace(tmp, _p(name))
            return
        except OSError:                                                 # noqa: PERF203
            if attempt == 4:
                raise
            time.sleep(0.5 + attempt)


def _read(name, default=None):
    try:
        with open(_p(name)) as fh:
            return json.load(fh)
    except Exception:                                                   # noqa: BLE001
        return default


def targets():
    """THE 60.  Frozen by `results/benchmark_manifest.json`; never filtered on anything."""
    return db.benchmark()


# ============================================================ stage: pool
_FRAG = {}


def fold_fragments(fold, n_folds=N_FOLDS):
    if fold not in _FRAG:
        _FRAG[fold] = list(dgm._fold_fragments(fold, n_folds))
    return _FRAG[fold]


def library_members(target_seq, fold, n_folds=N_FOLDS):
    """The leakage-safe library for one target: out-of-fold peptides + fold fragments.

    Identical rule to `s7.audit.build_pool_members`.  Folds assign whole IDENTITY CLUSTERS,
    so anything at or above `db.IDENTITY_THRESHOLD` to the target is in the target's own
    cluster, hence its own fold, hence excluded here; `stage_audit` verifies that by
    aligning every member rather than trusting the argument.
    """
    folds = db.folds(n_folds)
    peps = [q for q in db.load() if folds[q.seq] != fold and q.seq != target_seq]
    return peps, fold_fragments(fold, n_folds)


def windows_all(pool, n):
    """Every length-n window of every library member: CA, phi, psi, encoded seq, origin.

    Same iteration order as `s7.audit.windows_of` and `s8.generate._windows_all`, which is
    what makes the BLOSUM argsort -- and therefore the pool -- reproduce theirs exactly.
    """
    cas, phis, psis, seqs, src = [], [], [], [], []
    for q in pool:
        m = len(q.seq)
        if m < n:
            continue
        e = audit.encode(q.seq)
        for s in range(m - n + 1):
            cas.append(q.ca[s:s + n])
            phis.append(np.asarray(q.phi, float)[s:s + n])
            psis.append(np.asarray(q.psi, float)[s:s + n])
            seqs.append(e[s:s + n])
            src.append(q.pdb)
    return (np.stack(cas), np.stack(phis), np.stack(psis), np.stack(seqs),
            np.array(src, dtype=object))


_MODELS = {}
_GUARDED = [False]


def guard_esm():
    """Refuse to materialise the 1.5 GB `esm_cache.npz`.  Same guard as `s7.poolsize`.

    `dgm.train_fold` probes the feature width with ONE arbitrary TRAINING sequence, which
    is not in the 21 MB hot cache, and `esm_features.raw` would then load the full cache
    and take the box past its memory gate.  Every fold model is already on disk, so
    nothing trains and the probe's VALUES are never used -- only its shape, which fixes
    `d_in`.  So the probe is served zeros of the right shape and every real (cached)
    sequence goes through unchanged.  A BENCHMARK sequence missing from the hot cache is a
    hard error, never a silent zero.
    """
    if _GUARDED[0]:
        return
    import esm_features as ef
    for f in range(N_FOLDS):
        path = dgm._model_path(f, True, True, 0)
        if not os.path.exists(path):
            raise RuntimeError(f"missing fold model {path}; this stage must not train")
    hot = ef._load_small()
    real = {p.seq for p in targets()}
    _orig = ef.raw

    def raw(sequence):
        if sequence in hot:
            return _orig(sequence)
        if sequence in real:
            raise RuntimeError(f"benchmark sequence {sequence!r} absent from the ESM hot "
                               f"cache; refusing to fall back to the 1.5 GB cache")
        return (np.zeros((len(sequence), 1280), np.float32),
                np.zeros((len(sequence), len(sequence)), np.float32))

    ef.raw = raw
    _GUARDED[0] = True


def distogram_risk(seq, fold):
    """The fold model that never saw this target, applied to its sequence.

    Same call as `s7.debias.distogram_for`: `train_fold` LOADS the model off disk (all
    five exist) and trains nothing.  Only `_risk` is needed -- the shipped Bayes-risk
    lookup that `score_risk` gathers into.
    """
    guard_esm()
    if fold not in _MODELS:
        _MODELS[fold] = dgm.train_fold(fold, True, N_FOLDS, fragments=True, verbose=False)
    d = dgm.Distogram.for_target(seq, model=_MODELS[fold])
    return np.asarray(d._risk, float)


def _cpath(pdbid):
    return os.path.join(CACHE, f"{pdbid}.npz")


def build_pool(p, folds):
    """One benchmark target's complete record.  `rr`/`nat_ca`/`Dnat` are LABELS only."""
    fold = folds[p.seq]
    peps, frags = library_members(p.seq, fold)
    W, PHI, PSI, S, src = windows_all(peps + frags, p.n)
    sim = audit.B62[S, audit.encode(p.seq)[None, :]].sum(1)
    order = np.argsort(-sim, kind="stable")
    idx = order[:K]
    i, j = audit.pair_index(p.n)
    Wk = W[idx]
    D = audit.pair_dists(Wk, i, j)
    risk = distogram_risk(p.seq, fold)
    sc = debias.score_risk(risk, GRID, D)
    P = np.zeros((K, K), np.float32)
    for a in range(K):
        P[a] = audit.kabsch_rmsd_batch(Wk, Wk[a])
    return {"pdb": p.pdb, "n": int(p.n), "fold": int(fold), "seq": p.seq,
            "n_windows": int(len(W)), "n_peptides": len(peps), "n_fragments": len(frags),
            "W": Wk.astype(np.float32), "PHI": PHI[idx].astype(np.float32),
            "PSI": PSI[idx].astype(np.float32), "S": S[idx].astype(np.int8),
            "src": src[idx], "sim": sim[idx].astype(np.float32),
            "sc": sc.astype(np.float32), "D": D.astype(np.float32), "P": P,
            "rr": audit.kabsch_rmsd_batch(Wk, p.ca).astype(np.float32),
            "Dnat": audit.pair_dists(np.asarray(p.ca, float)[None], i, j)[0].astype(np.float32),
            "nat_ca": np.asarray(p.ca, np.float32)}


def load_pool(pdbid):
    z = np.load(_cpath(pdbid), allow_pickle=True)
    out = {"pdb": str(z["pdb"]), "n": int(z["n"]), "fold": int(z["fold"]),
           "seq": str(z["seq"]), "n_windows": int(z["n_windows"]),
           "n_peptides": int(z["n_peptides"]), "n_fragments": int(z["n_fragments"]),
           "src": z["src"], "S": z["S"].astype(int)}
    for k in ("W", "PHI", "PSI", "sim", "sc", "D", "P", "rr", "Dnat", "nat_ca"):
        out[k] = z[k].astype(float)
    return out


def deployable_view(c):
    """The ONLY object the synthesis and the relaxation are allowed to see.

    `rr`, `nat_ca` and `Dnat` are deliberately absent, so no downstream stage can read a
    native even by accident.  `stage_leak` NaN-poisons them in the cache and asserts every
    emitted structure is bit-identical, which is the test that this boundary is real.
    """
    return {k: c[k] for k in ("pdb", "n", "fold", "seq", "W", "PHI", "PSI", "S",
                              "sim", "sc", "P")}


def cached():
    if not os.path.isdir(CACHE):
        return []
    return sorted(f[:-4] for f in os.listdir(CACHE) if f.endswith(".npz"))


def stage_pool(verbose=True):
    os.makedirs(CACHE, exist_ok=True)
    if not start_gate():
        return
    folds = db.folds(N_FOLDS)
    tg = targets()
    todo = [p for p in tg if not os.path.exists(_cpath(p.pdb))]
    if verbose:
        print(f"{len(todo)}/{len(tg)} benchmark pools to build", flush=True)
    for k, p in enumerate(sorted(todo, key=lambda x: folds[x.seq])):
        if not gate(f"pool {p.pdb}"):
            return
        t0 = time.time()
        rec = build_pool(p, folds)
        tmp = _cpath(p.pdb) + ".tmp.npz"
        np.savez_compressed(tmp, **{kk: np.asarray(vv) for kk, vv in rec.items()})
        os.replace(tmp, _cpath(p.pdb))
        if verbose:
            print(f"[{k+1:3d}/{len(todo)}] {p.pdb:6} n={rec['n']:2d} fold={rec['fold']} "
                  f"windows={rec['n_windows']:6d} best@{K}={rec['rr'].min():.3f} "
                  f"shipped={float(rec['rr'][int(np.argmin(rec['sc']))]):.3f} "
                  f"({time.time()-t0:.0f}s) [mem {mem_pct()}%]", flush=True)
    if verbose:
        print(f"pools cached: {len(cached())}/{len(tg)}", flush=True)


# ============================================================ the window-identity control
#: Sub-window identity is NOT what the project's leakage filter controls.  `folds()`
#: assigns whole identity CLUSTERS and `_fold_fragments` aligns whole CHAINS, and
#: `peptide_db.identity` normalises by the LONGER sequence -- so a 60-residue fragment can
#: contain the target's exact 11-mer and still align at 0.52 to it, pass every filter, and
#: put an identity-1.0 window in the pool.  `stage_audit` measures that this happens on
#: 13 of the 60 benchmark targets, twice at identity 1.0.  This control rebuilds the pool
#: with every such window removed BEFORE the top-K cut -- so the pool is still K=500, made
#: of the next-best windows -- and runs the identical downstream.  It is a leakage control,
#: not an arm: the pre-registered headline is the uncontrolled one, and both are reported.
WIN_THR = 0.6                 #: the same identity threshold the rest of the project uses
CTRL_CACHE = os.path.join(HERE, "final_ctrl_cache")


def _ctrl_path(pdbid):
    return os.path.join(CTRL_CACHE, f"{pdbid}.npz")


def window_identity(target_seq, S_row):
    from peptide_db import identity
    return identity(target_seq, "".join(audit.AA[int(x)] for x in S_row))


def build_pool_clean(p, folds, thr=WIN_THR):
    """`build_pool` with every window at or above `thr` identity to the target removed.

    The exclusion happens BEFORE the top-K cut, walking down the BLOSUM order and skipping
    offending windows, so the control pool is still exactly K windows -- the next-best ones
    -- rather than a shortened version of the contaminated pool.  That is the pool the
    pipeline would have retrieved had the filter been applied at window level.
    """
    fold = folds[p.seq]
    peps, frags = library_members(p.seq, fold)
    W, PHI, PSI, S, src = windows_all(peps + frags, p.n)
    sim = audit.B62[S, audit.encode(p.seq)[None, :]].sum(1)
    order = np.argsort(-sim, kind="stable")
    keep, dropped = [], []
    for t in order:
        if window_identity(p.seq, S[t]) >= thr:
            dropped.append(int(t))
            continue
        keep.append(int(t))
        if len(keep) == K:
            break
    idx = np.asarray(keep, int)
    i, j = audit.pair_index(p.n)
    Wk = W[idx]
    D = audit.pair_dists(Wk, i, j)
    sc = debias.score_risk(distogram_risk(p.seq, fold), GRID, D)
    P = np.zeros((len(idx), len(idx)), np.float32)
    for a in range(len(idx)):
        P[a] = audit.kabsch_rmsd_batch(Wk, Wk[a])
    return {"pdb": p.pdb, "n": int(p.n), "fold": int(fold), "seq": p.seq,
            "n_windows": int(len(W)), "n_peptides": len(peps), "n_fragments": len(frags),
            "n_dropped_above_cut": len(dropped),
            "W": Wk.astype(np.float32), "PHI": PHI[idx].astype(np.float32),
            "PSI": PSI[idx].astype(np.float32), "S": S[idx].astype(np.int8),
            "src": src[idx], "sim": sim[idx].astype(np.float32),
            "sc": sc.astype(np.float32), "D": D.astype(np.float32), "P": P,
            "rr": audit.kabsch_rmsd_batch(Wk, p.ca).astype(np.float32),
            "Dnat": audit.pair_dists(np.asarray(p.ca, float)[None], i, j)[0].astype(np.float32),
            "nat_ca": np.asarray(p.ca, np.float32)}


def stage_ctrl(verbose=True):
    """The window-identity leakage control, end to end: pool -> synthesis -> AMBER."""
    import torsion_lib2 as tl2
    import amber_refine as ar
    os.makedirs(CTRL_CACHE, exist_ok=True)
    if not start_gate():
        return
    folds = db.folds(N_FOLDS)
    out = _read("final_ctrl.json") or {}
    tg = [p for p in targets() if p.pdb in set(cached())]
    for k, p in enumerate(tg):
        if p.pdb in out:
            continue
        if not gate(f"ctrl {p.pdb}"):
            return
        t0 = time.time()
        if os.path.exists(_ctrl_path(p.pdb)):
            z = np.load(_ctrl_path(p.pdb), allow_pickle=True)
            c = {"pdb": str(z["pdb"]), "n": int(z["n"]), "fold": int(z["fold"]),
                 "seq": str(z["seq"]), "S": z["S"].astype(int), "src": z["src"],
                 "n_dropped_above_cut": int(z["n_dropped_above_cut"])}
            for kk in ("W", "PHI", "PSI", "sim", "sc", "D", "P", "rr", "Dnat", "nat_ca"):
                c[kk] = z[kk].astype(float)
        else:
            c = build_pool_clean(p, folds)
            tmp = _ctrl_path(p.pdb) + ".tmp.npz"
            np.savez_compressed(tmp, **{kk: np.asarray(vv) for kk, vv in c.items()})
            os.replace(tmp, _ctrl_path(p.pdb))
        ca, phi, psi, fit, C, sub = synthesise(deployable_view(c))
        BB = geo.build_backbone_batch(np.asarray(phi)[None], np.asarray(psi)[None])
        cd = {kk: vv[0] for kk, vv in BB.items()}
        rep = tl2.PerResidueTorsion(c["seq"], tl2.library_for(c["seq"], 8, c["seq"]),
                                    chi_bits=False)
        full = None
        try:
            r = ar.refine_coords(c["seq"], rep, cd, k_restraint=AMBER_K,
                                 steps=AMBER_STEPS, components=True)
            full = np.asarray(r["ca"], float)
        except Exception as exc:                                        # noqa: BLE001
            print(f"  {p.pdb}: AMBER FAILED in control: {str(exc)[:160]}", flush=True)
        ar.clear_cache()
        # ---- natives enter ONLY here ----
        nat = c["nat_ca"]
        from s8.inband import sel_of
        rec = {"pdb": p.pdb, "n": c["n"], "fold": c["fold"],
               "n_dropped_above_cut": c["n_dropped_above_cut"],
               "synth": float(audit.kabsch_rmsd_batch(np.asarray(ca)[None], nat)[0]),
               "full": (float(audit.kabsch_rmsd_batch(full[None], nat)[0])
                        if full is not None else None),
               "shipped": sel_of(c["sc"], c["rr"]),
               "pool_best": float(c["rr"].min()),
               "wall": round(time.time() - t0, 1)}
        out[p.pdb] = rec
        _write("final_ctrl.json", out)
        if verbose:
            print(f"[{len(out):3d}/{len(tg)}] {p.pdb:6} dropped={rec['n_dropped_above_cut']:3d} "
                  f"synth={rec['synth']:.3f} full={rec['full']:.3f} "
                  f"pool={rec['pool_best']:.3f} ({rec['wall']:.0f}s) "
                  f"[mem {mem_pct()}%]", flush=True)
    if verbose and out:
        a = np.array([r["full"] for r in out.values() if r["full"] is not None])
        print(f"CONTROL  n={len(a)}  mean={a.mean():.4f}", flush=True)


# ============================================================ stage: synth (DEPLOYABLE)
def synthesise(v, m=M, pen_kind=PEN, lam=LAM):
    """Stages 2 and 3.  Nothing here reads a native -- `v` is a `deployable_view`.

    Returns (CA_arm, phi, psi, CA_fit, C_avg), where `fit` is the same projection with the
    prior switched off: the S8-11 incumbent, kept as the internal reference that the
    lambda = 0 rung of `lam_path` is asserted to reproduce exactly.
    """
    sub = np.argsort(np.asarray(v["sc"], float), kind="stable")[:m]
    Ps = np.asarray(v["P"], float)[np.ix_(sub, sub)]
    C = cc._avg_struct(np.asarray(v["W"], float)[sub], Ps)
    pen = pj.make_penalty(pen_kind, v["seq"], int(v["fold"]))
    path = pj.lam_path(C, pen, (0.0, lam), multi=True)
    fit, arm = path[0.0], path[lam]
    return arm[0], arm[1], arm[2], fit[0], C, sub


def stage_synth(verbose=True):
    if not start_gate():
        return
    out = _read("final_synth.json") or {}
    tg = cached()
    for k, pdbid in enumerate(tg):
        if pdbid in out:
            continue
        if not gate(f"synth {pdbid}"):
            return
        t0 = time.time()
        c = load_pool(pdbid)
        v = deployable_view(c)
        ca, phi, psi, fit, C, sub = synthesise(v)
        rec = {"pdb": pdbid, "n": c["n"], "fold": c["fold"], "seq": c["seq"],
               "ca": np.asarray(ca).tolist(), "phi": np.asarray(phi).tolist(),
               "psi": np.asarray(psi).tolist(), "fit_ca": np.asarray(fit).tolist(),
               "avg_ca": np.asarray(C).tolist(), "sub": np.asarray(sub, int).tolist(),
               "wall": round(time.time() - t0, 2)}
        # ---- natives enter ONLY here, after the structure is final ----
        nat = c["nat_ca"]
        rec["rmsd_arm"] = float(audit.kabsch_rmsd_batch(np.asarray(ca)[None], nat)[0])
        rec["rmsd_fit"] = float(audit.kabsch_rmsd_batch(np.asarray(fit)[None], nat)[0])
        rec["rmsd_avg"] = float(audit.kabsch_rmsd_batch(np.asarray(C)[None], nat)[0])
        out[pdbid] = rec
        _write("final_synth.json", out)
        if verbose:
            print(f"[{len(out):3d}/{len(tg)}] {pdbid:6} arm={rec['rmsd_arm']:.3f} "
                  f"fit={rec['rmsd_fit']:.3f} avg={rec['rmsd_avg']:.3f} "
                  f"({rec['wall']:.1f}s) [mem {mem_pct()}%]", flush=True)
    if verbose:
        a = np.array([r["rmsd_arm"] for r in out.values()])
        print(f"SYNTH  n={len(a)}  mean={a.mean():.4f}", flush=True)


# ============================================================ stage: amber (VALIDITY)
def stage_amber(verbose=True):
    """One restrained ff14SB/GBn2 minimisation per target.  ONE context, torn down after.

    S8-12 priced this exactly: k = 10-100 removes ~10^4 kcal/mol of builder strain for
    +0.011 to +0.012 A of CA-RMSD; S9-5 measured -0.0264 in its own harness.  It is a
    validity stage.  It is applied to the projected arm, which is what the pipeline emits.
    """
    import torsion_lib2 as tl2
    import amber_refine as ar
    if not start_gate():
        return
    syn = _read("final_synth.json") or {}
    out = _read("final_amber.json") or {}
    todo = [k for k in sorted(syn) if k not in out]
    if verbose:
        print(f"{len(todo)}/{len(syn)} relaxations to run", flush=True)
    for k, pdbid in enumerate(todo):
        if not gate(f"amber {pdbid}"):
            return
        r0 = syn[pdbid]
        seq = r0["seq"]
        phi = np.asarray(r0["phi"], float)
        psi = np.asarray(r0["psi"], float)
        BB = geo.build_backbone_batch(phi[None], psi[None])
        cd = {kk: vv[0] for kk, vv in BB.items()}
        tab = tl2.library_for(seq, 8, seq)             # holds the target's own seq out
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        t0 = time.time()
        try:
            r = ar.refine_coords(seq, rep, cd, k_restraint=AMBER_K, steps=AMBER_STEPS,
                                 components=True)
        except Exception as exc:                                        # noqa: BLE001
            out[pdbid] = {"pdb": pdbid, "err": str(exc)[:300]}
            _write("final_amber.json", out)
            print(f"  {pdbid}: AMBER FAILED: {str(exc)[:160]}", flush=True)
            ar.clear_cache()
            continue
        ca1 = np.asarray(r["ca"], float)
        bb = r.get("backbone", {})
        cm = r.get("components", {}) or {}
        rec = {"pdb": pdbid, "n": r0["n"], "seq": seq,
               "ca": ca1.tolist(),
               "e0": float(r["energy_initial"]), "e1": float(r["energy"]),
               "strain_after": float(cm.get("bond", 0.0) + cm.get("angle", 0.0)) if cm else None,
               "moved": float(r["restraint_rmsd"]),
               "wall": round(float(r["wall"]), 2)}
        if all(kk in bb for kk in ("N", "CA", "C")):
            ph1, ps1 = geo.extract_torsions(bb["N"], bb["CA"], bb["C"])
            rec["phi"] = np.asarray(ph1).tolist()
            rec["psi"] = np.asarray(ps1).tolist()
        # ---- natives enter ONLY here, after the physics is done ----
        c = load_pool(pdbid)
        rec["rmsd"] = float(audit.kabsch_rmsd_batch(ca1[None], c["nat_ca"])[0])
        rec["rmsd_before"] = r0["rmsd_arm"]
        rec["d_rmsd"] = rec["rmsd"] - r0["rmsd_arm"]
        ar.clear_cache()
        out[pdbid] = rec
        _write("final_amber.json", out)
        if verbose:
            print(f"[{k+1:3d}/{len(todo)}] {pdbid:6} {r0['rmsd_arm']:.3f} -> "
                  f"{rec['rmsd']:.3f}  dE={rec['e0']-rec['e1']:.3e}  "
                  f"moved={rec['moved']:.3f} ({time.time()-t0:.1f}s) "
                  f"[mem {mem_pct()}%]", flush=True)
    ok = [r for r in out.values() if "rmsd" in r]
    if verbose and ok:
        a = np.array([r["rmsd"] for r in ok])
        b = np.array([r["rmsd_before"] for r in ok])
        print(f"AMBER  n={len(a)}  {b.mean():.4f} -> {a.mean():.4f}  "
              f"d={a.mean()-b.mean():+.4f}", flush=True)


# ============================================================ statistics
def paired(a, b):
    """Paired difference a - b with a 95% CI, win/loss counts and the achieved SE."""
    d = np.asarray(a, float) - np.asarray(b, float)
    d = d[np.isfinite(d)]
    n = len(d)
    if n < 2:
        return None
    se = float(d.std(ddof=1) / math.sqrt(n))
    return {"mean_diff": float(d.mean()), "se": se,
            "ci95": [float(d.mean() - 1.96 * se), float(d.mean() + 1.96 * se)],
            "n_better": int((d < -1e-9).sum()), "n_worse": int((d > 1e-9).sum()),
            "n_tied": int((np.abs(d) <= 1e-9).sum()), "n": n}


def distribution(v):
    v = np.asarray(v, float)
    return {"n": int(len(v)), "mean": float(v.mean()), "median": float(np.median(v)),
            "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
            "min": float(v.min()), "max": float(v.max()),
            "frac_under_2.0": float((v < 2.0).mean()),
            "frac_under_1.5": float((v < 1.5).mean()),
            "frac_under_3.0": float((v < 3.0).mean())}


def concentration(d, drops=(5, 10, 20)):
    """How much of an aggregate gain is carried by its largest contributors.

    `d` is the per-target signed difference (negative = the system is better).  The share
    is the sum of the `k` most negative entries over the total, and the residual is the
    MEAN of `d` after removing them.  On the tuning instrument ten targets carried 61% of
    the incumbent's gain and it fell to -0.013 after dropping twenty; a mean carried by a
    handful of targets is much weaker than it looks, so this is reported unconditionally.
    """
    d = np.asarray(d, float)
    tot = float(d.sum())
    order = np.argsort(d)                       # most negative (largest gain) first
    out = {"mean": float(d.mean()), "total": tot, "n": int(len(d)),
           "sum_abs": float(np.abs(d).sum()), "n_gain": int((d < 0).sum()),
           "sum_gains": float(d[d < 0].sum()), "sum_losses": float(d[d > 0].sum())}
    #: A SHARE of a near-zero total is not a meaningful quantity: when the wins and the
    #: losses nearly cancel, dividing by their difference produces numbers in the hundreds
    #: of percent that say nothing.  The share is reported only when the total is at least
    #: a fifth of the total absolute movement; otherwise the honest statement is the
    #: absolute one -- how many Angstroms the top k carry, and what the mean becomes.
    out["share_is_meaningful"] = bool(abs(tot) >= 0.2 * out["sum_abs"])
    for k in drops:
        if k >= len(d):
            continue
        top = d[order[:k]]
        rest = d[order[k:]]
        out[f"top{k}_sum"] = float(top.sum())
        out[f"top{k}_share"] = (float(top.sum() / tot)
                                if out["share_is_meaningful"] and abs(tot) > 1e-12
                                else None)
        out[f"mean_ex_top{k}"] = float(rest.mean())
    return out


# ============================================================ stage: audit (LEAKAGE)
def stage_audit(verbose=True):
    """The leakage audit.  Counts and worst cases, never assurances.

    Four independent questions:
      1. Is every benchmark target scored by models and libraries that never saw it?
         -- fold discipline, checked per target rather than assumed.
      2. What is the WORST identity between a benchmark target and any library member
         (out-of-fold peptide or fold fragment) that its own pool was drawn from?
      3. What is the worst identity between a target and any POOL WINDOW's sequence?
         This is the sharper question: a window is a subsequence of a member, and identity
         normalises by the longer sequence, so a short window can score higher than its
         parent chain.  It is reported separately and in full.
      4. Does the benchmark intersect the 126-target tuning instrument or the 24-target
         dev set -- at pdb id, at sequence, and at identity CLUSTER?
    """
    from peptide_db import identity, max_possible_identity, _composition
    tg = targets()
    folds = db.folds(N_FOLDS)
    cl = db.clusters()
    tune = debias.tuning_targets()
    dev = db.dev_set(24)

    bset = {p.pdb for p in tg}
    out = {"n_benchmark": len(tg), "identity_threshold": db.IDENTITY_THRESHOLD,
           "intersections": {
               "benchmark_vs_tuning_pdb": sorted(bset & {p.pdb for p in tune}),
               "benchmark_vs_dev_pdb": sorted(bset & {p.pdb for p in dev}),
               "benchmark_vs_tuning_seq": sorted({p.seq for p in tg} & {p.seq for p in tune}),
               "benchmark_vs_dev_seq": sorted({p.seq for p in tg} & {p.seq for p in dev}),
               "benchmark_vs_tuning_cluster": sorted({cl[p.seq] for p in tg}
                                                     & {cl[p.seq] for p in tune}),
               "benchmark_vs_dev_cluster": sorted({cl[p.seq] for p in tg}
                                                  & {cl[p.seq] for p in dev}),
               "n_tuning": len(tune), "n_dev": len(dev)},
           "per_target": []}

    worst_mem = {"identity": -1.0}
    worst_win = {"identity": -1.0}
    n_mem_viol = n_win_viol = 0
    n_pairs = 0
    for k, p in enumerate(tg):
        fold = folds[p.seq]
        peps, frags = library_members(p.seq, fold)
        ct = _composition(p.seq)
        bm, bmw = 0.0, None
        for q in peps + frags:
            n_pairs += 1
            if max_possible_identity(p.seq, q.seq, ct) < bm:
                continue
            v = identity(p.seq, q.seq)
            if v > bm:
                bm, bmw = v, q.pdb
        if bm >= db.IDENTITY_THRESHOLD:
            n_mem_viol += 1
        if bm > worst_mem["identity"]:
            worst_mem = {"identity": float(bm), "target": p.pdb, "member": bmw}
        # every window that actually entered the K=500 pool
        bw, bww = 0.0, None
        if os.path.exists(_cpath(p.pdb)):
            c = load_pool(p.pdb)
            for b in range(len(c["S"])):
                w = "".join(audit.AA[int(x)] for x in c["S"][b])
                v = identity(p.seq, w)
                if v > bw:
                    bw, bww = v, str(c["src"][b])
            if bw >= db.IDENTITY_THRESHOLD:
                n_win_viol += 1
            if bw > worst_win["identity"]:
                worst_win = {"identity": float(bw), "target": p.pdb, "src": bww}
        # fold discipline, verified rather than assumed
        same_fold = [q.pdb for q in peps if folds[q.seq] == fold]
        out["per_target"].append({
            "pdb": p.pdb, "n": int(p.n), "fold": int(fold),
            "n_peptides": len(peps), "n_fragments": len(frags),
            "worst_member_identity": float(bm), "worst_member": bmw,
            "worst_window_identity": float(bw), "worst_window_src": bww,
            "self_in_library": p.seq in {q.seq for q in peps},
            "n_same_fold_in_library": len(same_fold)})
        if verbose and (k + 1) % 10 == 0:
            print(f"  audited {k+1}/{len(tg)}", flush=True)
    out["worst_member_identity"] = worst_mem
    out["worst_window_identity"] = worst_win
    out["n_member_violations"] = n_mem_viol
    out["n_window_violations"] = n_win_viol
    out["n_library_pairs_checked"] = n_pairs
    out["n_self_in_library"] = sum(r["self_in_library"] for r in out["per_target"])
    out["n_same_fold_in_library"] = sum(r["n_same_fold_in_library"]
                                        for r in out["per_target"])
    _write("final_audit.json", out)
    if verbose:
        _print_audit(out)
    return out


def _print_audit(out):
    ix = out["intersections"]
    print("\n=== LEAKAGE AUDIT ===")
    print(f"benchmark targets                      {out['n_benchmark']}")
    print(f"identity threshold                     {out['identity_threshold']}")
    print(f"library pairs aligned                  {out['n_library_pairs_checked']}")
    print(f"worst target-vs-LIBRARY-MEMBER identity {out['worst_member_identity']}")
    print(f"  violations (>= threshold)            {out['n_member_violations']}")
    print(f"worst target-vs-POOL-WINDOW identity    {out['worst_window_identity']}")
    print(f"  violations (>= threshold)            {out['n_window_violations']}")
    print(f"target's own sequence in its library    {out['n_self_in_library']}")
    print(f"same-fold peptides in any library       {out['n_same_fold_in_library']}")
    print(f"benchmark n {out['n_benchmark']} | tuning n {ix['n_tuning']} | dev n {ix['n_dev']}")
    for k in ("benchmark_vs_tuning_pdb", "benchmark_vs_dev_pdb",
              "benchmark_vs_tuning_seq", "benchmark_vs_dev_seq",
              "benchmark_vs_tuning_cluster", "benchmark_vs_dev_cluster"):
        print(f"  {k:34} {len(ix[k]):3d}  {ix[k][:6]}")


# ============================================================ stage: leak (NaN poison)
def stage_leak(n_targets=8, verbose=True):
    """Replace every native quantity with NaN and assert the pipeline is bit-identical.

    This is the test that the reporting boundary is real: `nat_ca`, `rr` and `Dnat` are
    poisoned in the loaded record, the WHOLE deployable path is re-run from the cached
    pool, and every emitted coordinate, torsion and filter index must match to 0.0.
    """
    tg = cached()[:n_targets]
    syn = _read("final_synth.json") or {}
    rows, worst = [], 0.0
    for pdbid in tg:
        c = load_pool(pdbid)
        c["nat_ca"] = np.full_like(c["nat_ca"], np.nan)
        c["rr"] = np.full_like(c["rr"], np.nan)
        c["Dnat"] = np.full_like(c["Dnat"], np.nan)
        v = deployable_view(c)
        ca, phi, psi, fit, C, sub = synthesise(v)
        ref = syn.get(pdbid)
        if ref is None:
            continue
        d = {"ca": float(np.abs(np.asarray(ca) - np.asarray(ref["ca"])).max()),
             "phi": float(np.abs(np.asarray(phi) - np.asarray(ref["phi"])).max()),
             "psi": float(np.abs(np.asarray(psi) - np.asarray(ref["psi"])).max()),
             "fit": float(np.abs(np.asarray(fit) - np.asarray(ref["fit_ca"])).max()),
             "avg": float(np.abs(np.asarray(C) - np.asarray(ref["avg_ca"])).max()),
             "sub": float(np.abs(np.asarray(sub, int)
                                 - np.asarray(ref["sub"], int)).max())}
        worst = max(worst, max(d.values()))
        rows.append({"pdb": pdbid, **d})
        if verbose:
            print(f"  {pdbid}: worst |diff| {max(d.values()):.3e}", flush=True)
        assert np.isfinite(np.asarray(ca)).all(), f"{pdbid}: poison leaked into the output"
    out = {"n_targets": len(rows), "worst_abs_diff": worst, "per_target": rows,
           "quantities": ["ca", "phi", "psi", "fit", "avg", "sub"]}
    _write("final_leak.json", out)
    if verbose:
        print(f"\nLEAK: {len(rows)} targets x 6 quantities, worst |diff| {worst:.3e}",
              flush=True)
    assert worst == 0.0, f"native quantities reach the deployable path: {worst}"
    return out


# ============================================================ stage: report
def _oracle_sel(c):
    """DIAGNOSTIC ONLY: rank the pool by L1 agreement with the NATIVE distance matrix."""
    s = np.abs(np.asarray(c["D"], float) - np.asarray(c["Dnat"], float)[None, :]).mean(1)
    from s8.inband import sel_of
    return sel_of(s, c["rr"])


def collect():
    """Every per-target number the report needs, in one pass over the caches."""
    from s8.inband import sel_of
    syn = _read("final_synth.json") or {}
    amb = _read("final_amber.json") or {}
    rows = []
    for pdbid in cached():
        if pdbid not in syn:
            continue
        c = load_pool(pdbid)
        s = syn[pdbid]
        a = amb.get(pdbid, {})
        sub = np.asarray(s["sub"], int)
        r = {"pdb": pdbid, "n": c["n"], "fold": c["fold"], "seq": c["seq"],
             "shipped": sel_of(c["sc"], c["rr"]),
             "pool_best": float(c["rr"].min()),
             "top75_best": float(c["rr"][sub].min()),
             "top75_mean": float(c["rr"][sub].mean()),
             "oracle_dist": _oracle_sel(c),
             "avg": s["rmsd_avg"], "fit": s["rmsd_fit"], "synth": s["rmsd_arm"],
             "full": a.get("rmsd"), "amber_moved": a.get("moved"),
             "amber_de": (a["e0"] - a["e1"]) if "e0" in a else None,
             "amber_strain_after": a.get("strain_after"),
             "phi": s["phi"], "psi": s["psi"], "ca": s["ca"],
             "ca_full": a.get("ca"), "phi_full": a.get("phi"), "psi_full": a.get("psi")}
        rows.append(r)
    return rows


def geometry_table(rows):
    """CA-CA bond, pseudo-angle, positive-phi and clash rate, on both emitted arms.

    The positive-phi column EXCLUDES `phi[0]` and `psi[n-1]`: `build_backbone_batch` never
    reads them, so at the end of a fit they hold whatever the optimiser left there (S9-2,
    finding 3).  `posphi_con_nongly` is the honest number; `posphi_all` is carried only for
    like-for-like comparison with `s8/audit8`'s table.
    """
    out = {}
    for tag, ck, pk, sk in (("synth", "ca", "phi", "psi"),
                            ("full", "ca_full", "phi_full", "psi_full")):
        g, t = [], []
        for r in rows:
            if r.get(ck) is None:
                continue
            g.append(pj.geometry_report(np.asarray(r[ck], float)))
            if r.get(pk) is not None:
                t.append(pj.torsion_report(r["seq"], np.asarray(r[pk], float),
                                           np.asarray(r[sk], float)))
        if not g:
            continue
        out[tag] = {"n": len(g), "geom": pj._agg(g)}
        out[tag]["geom_sd"] = {k: float(np.std([x[k] for x in g], ddof=1))
                               for k in g[0]}
        if t:
            out[tag]["tors"] = {k: float(np.nanmean([x[k] for x in t])) for k in t[0]}
            # pooled positive-phi rate over residues, not the mean of per-target rates
            npos = sum(x["n_pos_con_nongly"] for x in t)
            ntot = sum(x["n_con_nongly"] for x in t)
            out[tag]["posphi_con_nongly_pooled"] = float(npos / ntot) if ntot else None
            out[tag]["n_con_nongly_residues"] = int(ntot)
    # the natives and the library windows, as the reference columns
    gn, gw = [], []
    syn = _read("final_synth.json") or {}
    for r in rows:
        c = load_pool(r["pdb"])
        gn.append(pj.geometry_report(c["nat_ca"]))
        sub = np.asarray(syn[r["pdb"]]["sub"], int)
        Ps = c["P"][np.ix_(sub, sub)]
        gw.append(pj.geometry_report(c["W"][sub][int(np.argmin(cc._medoid(Ps)))]))
    out["native"] = {"n": len(gn), "geom": pj._agg(gn),
                     "geom_sd": {k: float(np.std([x[k] for x in gn], ddof=1)) for k in gn[0]}}
    out["window"] = {"n": len(gw), "geom": pj._agg(gw),
                     "geom_sd": {k: float(np.std([x[k] for x in gw], ddof=1)) for k in gw[0]}}
    return out


def _replication(arms, n):
    """Did the tuning result replicate -- as a MEAN, and as the EFFECT it stands on?

    Two separate questions, and they have different answers, so they are separated here.
    The absolute mean measures how hard the target set is.  The paired gain against the
    shipped distogram argmin measures the architecture.  A benchmark whose absolute mean
    is lower but whose paired gain is zero has NOT replicated the result, and saying so is
    the whole point of a held-out pass.
    """
    p = paired(arms["full"], arms["shipped"])
    ps = paired(arms["synth"], arms["shipped"])
    lo, hi = p["ci95"]
    return {
        "tuning_mean": TUNING_MEAN, "dev_mean": DEV_MEAN,
        "tuning_baseline": TUNING_BASELINE, "tuning_pool_best": TUNING_POOL_BEST,
        "tuning_gain": TUNING_GAIN, "tuning_gain_ci": list(TUNING_GAIN_CI),
        "dev_gain": DEV_GAIN,
        "benchmark_full": float(arms["full"].mean()),
        "benchmark_synth": float(arms["synth"].mean()),
        "benchmark_baseline": float(arms["shipped"].mean()),
        "benchmark_pool_best": float(arms["pool_best"].mean()),
        "d_vs_tuning": float(arms["full"].mean() - TUNING_MEAN),
        "d_vs_tuning_ablation": float(arms["synth"].mean() - TUNING_MEAN),
        "d_vs_dev": float(arms["full"].mean() - DEV_MEAN),
        "se": float(arms["full"].std(ddof=1) / math.sqrt(n)),
        "benchmark_gain_full": p["mean_diff"], "benchmark_gain_full_ci": p["ci95"],
        "benchmark_gain_synth": ps["mean_diff"], "benchmark_gain_synth_ci": ps["ci95"],
        "paired_se": p["se"],
        #: the decisive test: could this pass have SEEN the tuning effect had it been
        #: there?  n * SE gives the resolution; the CI containing or excluding
        #: TUNING_GAIN gives the verdict.
        "tuning_gain_in_benchmark_ci": bool(lo <= TUNING_GAIN <= hi),
        "tuning_gain_in_benchmark_ci_synth": bool(ps["ci95"][0] <= TUNING_GAIN
                                                  <= ps["ci95"][1]),
        "gain_resolvable_at_2se": bool(abs(TUNING_GAIN) >= 2.0 * p["se"]),
        "tuning_gain_in_se_units": float(abs(TUNING_GAIN) / p["se"]),
        "verdict": ("the mean replicates and the EFFECT does not"
                    if not (lo <= TUNING_GAIN <= hi) else "the effect replicates")}


def _control_block():
    """The window-identity leakage control, folded into the report if it has been run."""
    ct = _read("final_ctrl.json")
    syn = _read("final_synth.json") or {}
    amb = _read("final_amber.json") or {}
    if not ct:
        return None
    au = _read("final_audit.json") or {}
    hot = {r["pdb"] for r in au.get("per_target", [])
           if r.get("worst_window_identity", 0.0) >= WIN_THR}
    keys = [k for k in sorted(ct)
            if ct[k].get("full") is not None and k in amb and "rmsd" in amb[k]]
    if not keys:
        return None
    cf = np.array([ct[k]["full"] for k in keys])
    mf = np.array([amb[k]["rmsd"] for k in keys])
    cs = np.array([ct[k]["synth"] for k in keys])
    ms = np.array([syn[k]["rmsd_arm"] for k in keys])
    ch = np.array([ct[k]["shipped"] for k in keys])
    out = {"n": len(keys), "threshold": WIN_THR,
           "n_targets_with_a_contaminated_window": len(hot),
           "n_windows_dropped_total": int(sum(ct[k]["n_dropped_above_cut"]
                                              for k in keys)),
           "control_full_mean": float(cf.mean()), "main_full_mean": float(mf.mean()),
           "control_synth_mean": float(cs.mean()), "main_synth_mean": float(ms.mean()),
           "control_shipped_mean": float(ch.mean()),
           "full_control_vs_main": paired(cf, mf),
           "synth_control_vs_main": paired(cs, ms),
           "control_full_vs_control_shipped": paired(cf, ch)}
    m = np.array([k in hot for k in keys])
    if m.any():
        out["on_contaminated_targets"] = {
            "n": int(m.sum()),
            "control_full": float(cf[m].mean()), "main_full": float(mf[m].mean()),
            "paired": paired(cf[m], mf[m]),
            "per_target": [{"pdb": keys[i], "main": float(mf[i]),
                            "control": float(cf[i]),
                            "dropped": ct[keys[i]]["n_dropped_above_cut"]}
                           for i in range(len(keys)) if m[i]]}
    if (~m).any():
        out["on_clean_targets"] = {"n": int((~m).sum()),
                                   "paired": paired(cf[~m], mf[~m])}
    return out


def stage_report(verbose=True):
    rows = collect()
    done = [r for r in rows if r["full"] is not None]
    arms = {"full": np.array([r["full"] for r in done], float),
            "synth": np.array([r["synth"] for r in done], float),
            "fit": np.array([r["fit"] for r in done], float),
            "avg": np.array([r["avg"] for r in done], float),
            "shipped": np.array([r["shipped"] for r in done], float),
            "pool_best": np.array([r["pool_best"] for r in done], float),
            "top75_best": np.array([r["top75_best"] for r in done], float),
            "oracle_dist": np.array([r["oracle_dist"] for r in done], float)}
    out = {"n": len(done), "n_pools": len(cached()),
           "preregistration": {"K": K, "m": M, "penalty": PEN, "lam": LAM,
                               "amber_k": AMBER_K, "amber_steps": AMBER_STEPS,
                               "multi_start": True, "n_folds": N_FOLDS,
                               "tuning_mean": TUNING_MEAN, "dev_mean": DEV_MEAN},
           "dist": {k: distribution(v) for k, v in arms.items()},
           "paired": {
               "full_vs_shipped": paired(arms["full"], arms["shipped"]),
               "synth_vs_shipped": paired(arms["synth"], arms["shipped"]),
               "full_vs_synth": paired(arms["full"], arms["synth"]),
               "synth_vs_fit": paired(arms["synth"], arms["fit"]),
               "full_vs_poolbest": paired(arms["full"], arms["pool_best"]),
               "shipped_vs_poolbest": paired(arms["shipped"], arms["pool_best"]),
               "oracle_vs_shipped": paired(arms["oracle_dist"], arms["shipped"])},
           "concentration": {
               "full_vs_shipped": concentration(arms["full"] - arms["shipped"]),
               "synth_vs_shipped": concentration(arms["synth"] - arms["shipped"])},
           "selection_gap": {
               "full_minus_pool_best": float(arms["full"].mean() - arms["pool_best"].mean()),
               "full_minus_top75_best": float(arms["full"].mean() - arms["top75_best"].mean()),
               "shipped_minus_pool_best": float(arms["shipped"].mean()
                                                - arms["pool_best"].mean())},
           "replication": _replication(arms, len(done)),
           "geometry": geometry_table(done),
           "window_identity_control": _control_block(),
           "leak": _read("final_leak.json", {}).get("worst_abs_diff"),
           "audit_summary": {k: (_read("final_audit.json") or {}).get(k) for k in
                             ("worst_member_identity", "n_member_violations",
                              "worst_window_identity", "n_window_violations",
                              "n_library_pairs_checked", "n_self_in_library",
                              "n_same_fold_in_library", "intersections")},
           "amber": {
               "mean_energy_drop": float(np.mean([r["amber_de"] for r in done
                                                  if r["amber_de"] is not None])),
               "median_energy_drop": float(np.median([r["amber_de"] for r in done
                                                      if r["amber_de"] is not None])),
               "mean_moved": float(np.mean([r["amber_moved"] for r in done
                                            if r["amber_moved"] is not None])),
               "mean_strain_after": float(np.mean([r["amber_strain_after"] for r in done
                                                   if r["amber_strain_after"] is not None]))},
           "by_length": {},
           "per_target": [{k: r[k] for k in
                           ("pdb", "n", "fold", "seq", "shipped", "pool_best",
                            "top75_best", "oracle_dist", "avg", "fit", "synth", "full",
                            "amber_moved")} for r in done]}
    for n in sorted({r["n"] for r in done}):
        m = np.array([r["n"] == n for r in done])
        out["by_length"][int(n)] = {"n_targets": int(m.sum()),
                                    "full": float(arms["full"][m].mean()),
                                    "shipped": float(arms["shipped"][m].mean()),
                                    "pool_best": float(arms["pool_best"][m].mean())}
    out["by_fold"] = {}
    for f in sorted({r["fold"] for r in done}):
        m = np.array([r["fold"] == f for r in done])
        out["by_fold"][int(f)] = {"n_targets": int(m.sum()),
                                  "full": float(arms["full"][m].mean()),
                                  "shipped": float(arms["shipped"][m].mean()),
                                  "d": float((arms["full"] - arms["shipped"])[m].mean())}
    _write("final_report.json", out)
    if verbose:
        _print_report(out)
    return out


def _print_report(o):
    p = o["preregistration"]
    print("\n" + "=" * 78)
    print("SPRINT 9 FINAL EVALUATION -- 60-TARGET HELD-OUT BENCHMARK, ONE PASS")
    print("=" * 78)
    print(f"PRE-REGISTERED: K={p['K']} pool -> distogram top-{p['m']} -> coordinate "
          f"average -> manifold projection")
    print(f"                (multi-start, {p['penalty']}@{p['lam']}) -> AMBER k="
          f"{p['amber_k']:g}, steps={p['amber_steps']} (converged)")
    print(f"n = {o['n']} targets\n")
    print(f"{'arm':22} {'mean':>8} {'median':>8} {'sd':>7} {'min':>7} {'max':>7} "
          f"{'<2.0':>6} {'<1.5':>6}")
    order = [("full", "FULL SYSTEM (+AMBER)"), ("synth", "stage 3 only (ablation)"),
             ("fit", "  no torsion prior"), ("avg", "  raw average (illegal)"),
             ("shipped", "shipped distogram argmin"), ("top75_best", "top-75 best (ceil)"),
             ("pool_best", "pool best (ceiling)"), ("oracle_dist", "distance oracle (diag)")]
    for k, name in order:
        d = o["dist"][k]
        print(f"{name:22} {d['mean']:8.4f} {d['median']:8.4f} {d['sd']:7.4f} "
              f"{d['min']:7.3f} {d['max']:7.3f} {d['frac_under_2.0']:6.3f} "
              f"{d['frac_under_1.5']:6.3f}")
    print("\n--- paired comparisons (negative = the first arm is better) ---")
    for k, v in o["paired"].items():
        if v is None:
            continue
        print(f"{k:24} d={v['mean_diff']:+.4f}  95% CI [{v['ci95'][0]:+.4f}, "
              f"{v['ci95'][1]:+.4f}]  SE {v['se']:.4f}  W/L {v['n_better']}/{v['n_worse']}")
    print("\n--- concentration of the gain over the shipped pipeline ---")
    for k, c in o["concentration"].items():
        print(f"{k}: mean {c['mean']:+.4f}  total {c['total']:+.3f} A over {c['n']} "
              f"targets; gains {c['sum_gains']:+.3f} on {c['n_gain']}, losses "
              f"{c['sum_losses']:+.3f} on {c['n'] - c['n_gain']}")
        if not c["share_is_meaningful"]:
            print("   [share of the total is NOT reported: the wins and losses nearly "
                  "cancel, so a percentage of their difference is not a real quantity]")
        for j in (5, 10, 20):
            if f"top{j}_sum" not in c:
                continue
            sh = (f"{c[f'top{j}_share']*100:6.1f}% of the total"
                  if c[f"top{j}_share"] is not None else
                  f"{c[f'top{j}_sum']:+.3f} A")
            print(f"   top-{j:2d} carry {sh}; mean after dropping them "
                  f"{c[f'mean_ex_top{j}']:+.4f}")
    print("\n--- selection gap ---")
    for k, v in o["selection_gap"].items():
        print(f"{k:28} {v:+.4f}")
    print("\n--- geometry validity ---")
    g = o["geometry"]
    ks = ["step_mean", "step_sd", "pseudoangle_mean", "frac_pseudoangle_out_of_75_150",
          "min_nonlocal_CA", "frac_nonlocal_under_4A", "rg"]
    cols = [c for c in ("full", "synth", "window", "native") if c in g]
    print(f"{'quantity':34} " + "".join(f"{c:>12}" for c in cols))
    for k in ks:
        print(f"{k:34} " + "".join(f"{g[c]['geom'][k]:12.4f}" for c in cols))
    print(f"{'CA-CA bond SD across targets':34} "
          + "".join(f"{g[c]['geom_sd']['step_mean']:12.2e}" for c in cols))
    for c in cols:
        if "tors" in g.get(c, {}):
            t = g[c]["tors"]
            print(f"  {c}: pos-phi non-Gly excl. phi[0]/psi[n-1] = "
                  f"{t['posphi_con_nongly']:.4f} (per-target mean), "
                  f"{g[c]['posphi_con_nongly_pooled']:.4f} pooled over "
                  f"{g[c]['n_con_nongly_residues']} residues; all-residue "
                  f"{t['posphi_all']:.4f}")
    a = o["amber"]
    print(f"\n--- AMBER validity stage ---")
    print(f"mean energy removed  {a['mean_energy_drop']:.3e} kcal/mol   "
          f"median {a['median_energy_drop']:.3e}")
    print(f"mean CA/N/C moved    {a['mean_moved']:.3f} A   "
          f"mean residual bond+angle strain {a['mean_strain_after']:.2f} kcal/mol")
    c = o.get("window_identity_control")
    if c:
        print(f"\n--- window-identity leakage control (drop pool windows >= "
              f"{c['threshold']} identity to the target, re-retrieve to K={K}) ---")
        print(f"{c['n_targets_with_a_contaminated_window']} of {c['n']} targets had at "
              f"least one such window; {c['n_windows_dropped_total']} windows dropped")
        print(f"full system   main {c['main_full_mean']:.4f} -> control "
              f"{c['control_full_mean']:.4f}   d={c['full_control_vs_main']['mean_diff']:+.4f} "
              f"[{c['full_control_vs_main']['ci95'][0]:+.4f}, "
              f"{c['full_control_vs_main']['ci95'][1]:+.4f}]")
        print(f"stage 3 only  main {c['main_synth_mean']:.4f} -> control "
              f"{c['control_synth_mean']:.4f}   d={c['synth_control_vs_main']['mean_diff']:+.4f} "
              f"[{c['synth_control_vs_main']['ci95'][0]:+.4f}, "
              f"{c['synth_control_vs_main']['ci95'][1]:+.4f}]")
        if "on_contaminated_targets" in c:
            h = c["on_contaminated_targets"]
            print(f"on the {h['n']} contaminated targets: {h['main_full']:.4f} -> "
                  f"{h['control_full']:.4f}  d={h['paired']['mean_diff']:+.4f} "
                  f"[{h['paired']['ci95'][0]:+.4f}, {h['paired']['ci95'][1]:+.4f}]")
            for t in h["per_target"]:
                print(f"    {t['pdb']:6} dropped={t['dropped']:3d} "
                      f"{t['main']:.3f} -> {t['control']:.3f}")
        if "on_clean_targets" in c:
            k = c["on_clean_targets"]["paired"]
            print(f"on the {c['on_clean_targets']['n']} clean targets (must be exactly "
                  f"zero): d={k['mean_diff']:+.6f}")
    r = o["replication"]
    print(f"\n--- replication ---")
    print(f"{'':22}{'tuning n=126':>14}{'dev n=24':>12}{'BENCHMARK n=' + str(o['n']):>16}")
    print(f"{'system mean':22}{r['tuning_mean']:14.3f}{r['dev_mean']:12.3f}"
          f"{r['benchmark_full']:16.4f}")
    print(f"{'shipped baseline':22}{r['tuning_baseline']:14.3f}{'--':>12}"
          f"{r['benchmark_baseline']:16.4f}")
    print(f"{'pool best (ceiling)':22}{r['tuning_pool_best']:14.3f}{'--':>12}"
          f"{r['benchmark_pool_best']:16.4f}")
    print(f"{'PAIRED GAIN vs shipped':22}{r['tuning_gain']:+14.3f}{r['dev_gain']:+12.3f}"
          f"{r['benchmark_gain_full']:+16.4f}")
    print(f"benchmark - tuning mean = {r['d_vs_tuning']:+.4f} (full) / "
          f"{r['d_vs_tuning_ablation']:+.4f} (ablation);  mean SE {r['se']:.4f}")
    print(f"benchmark paired gain 95% CI [{r['benchmark_gain_full_ci'][0]:+.4f}, "
          f"{r['benchmark_gain_full_ci'][1]:+.4f}], paired SE {r['paired_se']:.4f}")
    print(f"the tuning gain {r['tuning_gain']:+.3f} is {r['tuning_gain_in_se_units']:.1f} "
          f"paired SE from zero here, and is "
          f"{'INSIDE' if r['tuning_gain_in_benchmark_ci'] else 'OUTSIDE'} the "
          f"benchmark's CI")
    print(f"VERDICT: {r['verdict']}")
    print("\n--- per target ---")
    print(f"{'pdb':6} {'n':>3} {'fold':>4} {'shipped':>8} {'synth':>8} {'FULL':>8} "
          f"{'pool':>7} {'oracle':>7}")
    for t in sorted(o["per_target"], key=lambda x: x["full"]):
        print(f"{t['pdb']:6} {t['n']:3d} {t['fold']:4d} {t['shipped']:8.3f} "
              f"{t['synth']:8.3f} {t['full']:8.3f} {t['pool_best']:7.3f} "
              f"{t['oracle_dist']:7.3f}")


# ============================================================ main
STAGES = {"pool": stage_pool, "synth": stage_synth, "amber": stage_amber,
          "audit": stage_audit, "ctrl": stage_ctrl, "leak": stage_leak,
          "report": stage_report}


def main(argv):
    if not argv or argv[0] not in STAGES:
        print(__doc__)
        print(f"stages: {' '.join(STAGES)}")
        return 1
    STAGES[argv[0]]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
