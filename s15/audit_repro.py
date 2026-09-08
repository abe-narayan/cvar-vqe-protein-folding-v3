"""S15 AUDIT / Part A -- reproduce the instrument from primary inputs.

Stages, each writes its own JSON so a partial run is never mistaken for a complete one:

    python -m s15.audit_repro db      # rescan pdbs/ + pdbs_ext/ -> compare to peptide_db.npz
    python -m s15.audit_repro univ    # rebuild s8/generate_univ/<pdb>.npz for N targets
    python -m s15.audit_repro disto   # recompute s12/cache/disto_<pdb>.npz from the model
    python -m s15.audit_repro deps    # which of the five constants touches which artefact

Every comparison is bit-level where the stored dtype allows it, otherwise max-abs-diff
against a stated tolerance.
"""
from __future__ import annotations
import os, sys, json, time, glob, hashlib
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s15", "results")
os.makedirs(OUT, exist_ok=True)


def _w(name, obj):
    p = os.path.join(OUT, name if name.endswith(".json") else name + ".json")
    with open(p, "w") as fh:
        json.dump(obj, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print("wrote", p, flush=True)
    return p


def sha(path, nbytes=None):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(1 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()[:16]


def _cmp(a, b):
    """Return (bit_identical, max_abs_diff, n_diff) for two arrays."""
    a = np.asarray(a); b = np.asarray(b)
    if a.shape != b.shape:
        return False, float("inf"), -1
    if a.dtype.kind in "SU" or b.dtype.kind in "SU" or a.dtype == object:
        eq = bool(np.all(a == b))
        return eq, 0.0 if eq else float("inf"), int((a != b).sum())
    bit = bool(np.array_equal(a, b))
    if a.dtype.kind == "b":
        return bit, 0.0 if bit else 1.0, int((a != b).sum())
    d = np.abs(a.astype(np.float64) - b.astype(np.float64))
    return bit, float(d.max()) if d.size else 0.0, int((a != b).sum())


# ------------------------------------------------------------------ stage: db
def stage_db():
    """Rescan the primary PDB files and compare to the committed peptide_db.npz."""
    import peptide_db as db
    t0 = time.time()
    recs = db._scan()
    dt = time.time() - t0
    z = np.load(db.CACHE, allow_pickle=True)
    old = list(z["records"])
    res = {"n_scanned": len(recs), "n_cached": len(old), "scan_seconds": round(dt, 1),
           "cache_sha256_16": sha(db.CACHE),
           "dirs": [os.path.basename(d) for d in db.DIRS],
           "n_pdb_files": {os.path.basename(d): len(glob.glob(os.path.join(d, "*.pdb"))) for d in db.DIRS}}
    bynew = {r["seq"]: r for r in recs}
    byold = {r["seq"]: r for r in old}
    res["seq_sets_equal"] = sorted(bynew) == sorted(byold)
    res["only_new"] = sorted(set(bynew) - set(byold))[:20]
    res["only_cached"] = sorted(set(byold) - set(bynew))[:20]
    worst = {"ca": 0.0, "phi": 0.0, "psi": 0.0}
    nbit = {"ca": 0, "phi": 0, "psi": 0}
    npdb_diff = []
    common = sorted(set(bynew) & set(byold))
    for s in common:
        for k in ("ca", "phi", "psi"):
            bit, mx, _ = _cmp(bynew[s][k], byold[s][k])
            worst[k] = max(worst[k], mx)
            nbit[k] += int(bit)
        if bynew[s]["pdb"] != byold[s]["pdb"]:
            npdb_diff.append((byold[s]["pdb"], bynew[s]["pdb"]))
    res["n_common"] = len(common)
    res["max_abs_diff"] = worst
    res["n_bit_identical"] = nbit
    res["pdb_id_disagreements"] = npdb_diff[:20]
    res["verdict"] = ("BIT-IDENTICAL" if (res["seq_sets_equal"] and not npdb_diff
                                          and all(nbit[k] == len(common) for k in nbit))
                      else "DRIFT")
    return _w("audit_repro_db", res)


# ---------------------------------------------------------------- stage: univ
def stage_univ(pdbs=None, nmax=3):
    """Rebuild the window universe for a few targets, from the library, and compare."""
    import math
    import distogram as dgm
    import peptide_db as db
    from s7 import audit, debias
    from s8 import generate as g8

    tg = {p.pdb: p for p in debias.tuning_targets()}
    if pdbs is None:
        pdbs = sorted(tg)[:nmax]
    folds = db.folds(5)
    memo = {}
    _orig = dgm._fold_fragments

    def frags(fold, n_folds=5, threshold=db.IDENTITY_THRESHOLD):
        key = (fold, n_folds, threshold)
        if key not in memo:
            memo[key] = _orig(fold, n_folds, threshold)
        return memo[key]

    dgm._fold_fragments = frags
    rows = []
    try:
        for pdb in pdbs:
            p = tg[pdb]
            t0 = time.time()
            fold = folds[p.seq]
            peps = [q for q in db.load() if folds[q.seq] != fold and q.seq != p.seq]
            fr = list(frags(fold, 5))
            W, PH, PS, S, org = g8._windows_all(peps + fr, p.n,
                                                [True] * len(peps) + [False] * len(fr))
            sim = audit.B62[S, audit.encode(p.seq)[None, :]].sum(1)
            order = np.argsort(-sim, kind="stable")
            new = {"W": W.astype(np.float32), "PHI": PH.astype(np.float16),
                   "PSI": PS.astype(np.float16), "S": S.astype(np.int8), "org": org,
                   "sim": sim.astype(np.float32), "order": order.astype(np.int32),
                   "rr": audit.kabsch_rmsd_batch(W, p.ca).astype(np.float32),
                   "nat_ca": np.asarray(p.ca, np.float32)}
            path = os.path.join(ROOT, "s8", "generate_univ", f"{pdb}.npz")
            z = np.load(path, allow_pickle=True)
            row = {"pdb": pdb, "n": int(p.n), "fold": int(fold),
                   "n_windows_new": int(len(W)), "n_windows_cached": int(len(z["W"])),
                   "seconds": round(time.time() - t0, 1), "sha256_16": sha(path)}
            for k in new:
                bit, mx, nd = _cmp(new[k], z[k])
                row[k] = {"bit": bit, "max_abs_diff": mx, "n_elem_diff": nd}
            # the derived constants for this target
            K, M, BAND = 500, 75, 1.5
            o = np.asarray(z["order"], int)[:K]
            row["pool_best_cached"] = float(np.asarray(z["rr"], np.float64)[o].min())
            o2 = np.asarray(new["order"], int)[:K]
            row["pool_best_rebuilt"] = float(np.asarray(new["rr"], np.float64)[o2].min())
            row["pool_best_delta"] = row["pool_best_rebuilt"] - row["pool_best_cached"]
            rows.append(row)
            print(json.dumps(row), flush=True)
            _w("audit_repro_univ", {"rows": rows, "n_expected": len(pdbs),
                                    "complete": len(rows) == len(pdbs)})
    finally:
        dgm._fold_fragments = _orig
    return _w("audit_repro_univ", {"rows": rows, "n_expected": len(pdbs),
                                   "complete": len(rows) == len(pdbs)})


# --------------------------------------------------------------- stage: disto
def stage_disto(nmax=3, pdbs=None):
    """Recompute the leave-fold-out distogram from the model and compare to s12/cache."""
    from s12 import instrument as I
    from core import pipeline as pl
    from core import predict as dgm
    tg = I.targets()
    if pdbs is None:
        pdbs = [t["pdb"] for t in tg[:nmax]]
    byp = {t["pdb"]: t for t in tg}
    rows = []
    for pdb in pdbs:
        t = byp[pdb]
        path = os.path.join(I.CACHE, f"disto_{pdb}.npz")
        z = np.load(path)
        t0 = time.time()
        pl.guard_esm([t["seq"]])
        model = pl.fold_model(int(t["fold"]))
        d = dgm.Distogram.for_target(t["seq"], model=model)
        new = {"prob": d.prob.astype(np.float32), "i": d.i, "j": d.j,
               "expected": d.expected.astype(np.float32), "sd": d.sd.astype(np.float32),
               "risk": np.asarray(d._risk, np.float32)}
        row = {"pdb": pdb, "fold": int(t["fold"]), "seconds": round(time.time() - t0, 1),
               "sha256_16": sha(path)}
        for k in new:
            bit, mx, nd = _cmp(new[k], z[k])
            row[k] = {"bit": bit, "max_abs_diff": mx, "n_elem_diff": nd}
        # does the shipped argmin move?
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        i, j = I.pair_index(t["n"])
        D = I.pair_dists(u["W"][p], i, j).astype(np.float32).astype(float)
        old_sc = I.shipped_score({k: z[k] for k in z.files}, D)
        newd = dict({k: z[k] for k in z.files}); newd.update(new)
        new_sc = I.shipped_score(newd, D)
        rr = u["rr"][p]
        row["shipped_rmsd_cached"] = float(rr[int(np.argmin(old_sc))])
        row["shipped_rmsd_recomputed"] = float(rr[int(np.argmin(new_sc))])
        row["argmin_same"] = int(np.argmin(old_sc)) == int(np.argmin(new_sc))
        rows.append(row)
        print(json.dumps(row), flush=True)
        _w("audit_repro_disto", {"rows": rows, "n_expected": len(pdbs),
                                 "complete": len(rows) == len(pdbs)})
    return _w("audit_repro_disto", {"rows": rows, "n_expected": len(pdbs),
                                    "complete": len(rows) == len(pdbs)})


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "db"
    arg = sys.argv[2] if len(sys.argv) > 2 else None
    if cmd == "db":
        stage_db()
    elif cmd == "univ":
        if arg == "all":
            from s7 import debias
            stage_univ(pdbs=sorted(p.pdb for p in debias.tuning_targets()))
        else:
            stage_univ(pdbs=arg.split(",") if arg else None)
    elif cmd == "disto":
        stage_disto(pdbs=arg.split(",") if arg else None)
    else:
        raise SystemExit("unknown stage " + cmd)
