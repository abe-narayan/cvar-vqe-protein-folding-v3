"""FAIL18 forensics, step 4: NMR ensemble analysis (label-convention audit).

For every tuning target, parse ALL deposited models and measure
  - ensemble spread (mean pairwise CA-RMSD across deposited models)
  - RMSD of model-1 (the project's declared endpoint) to the ensemble medoid
  - the emitted structure's RMSD to model-1 (declared metric, UNCHANGED) and to the
    BEST deposited model (secondary diagnostic only)
  - the pool best against model-1 vs against any model

This is a DIAGNOSTIC arm.  The project's declared endpoint (model-1 CA-RMSD) is
reported first everywhere and is not replaced.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
from s12 import instrument as I
import protein_geometry as geo

PDBDIRS = [os.path.join(ROOT, "pdbs_ext"), os.path.join(ROOT, "pdbs")]


def path_for(pdb):
    for d in PDBDIRS:
        p = os.path.join(d, f"{pdb}.pdb")
        if os.path.exists(p):
            return p
    return None


def ensemble_ca(pdb):
    p = path_for(pdb)
    models = geo.parse_pdb_ensemble(p)
    return [np.asarray(m[2], float) for m in models]     # CA of each model


def main():
    tg = I.targets()
    out = []
    for k, t in enumerate(tg):
        pdb, n = t["pdb"], t["n"]
        try:
            E = ensemble_ca(pdb)
        except Exception as ex:
            out.append(dict(pdb=pdb, n=n, error=str(ex))); continue
        E = [e for e in E if len(e) == n]
        if not E:
            out.append(dict(pdb=pdb, n=n, error="length mismatch")); continue
        M = np.stack(E)                                   # (m, n, 3)
        m = len(M)
        u = I.load_univ(pdb); p = I.pool_idx(u)
        nat1 = u["nat_ca"]
        rec = I.shipped_record(pdb)
        emit = np.asarray(rec.get("amber_ca", rec["ca"]), float)

        # sanity: model 1 in the npz must be model 1 on disk
        d_check = I.ca_rmsd(M[0], nat1)

        if m > 1:
            P = I.pairwise_rmsd(M)
            spread = float(P[np.triu_indices(m, 1)].mean())
            med = int(np.argmin(P.mean(1)))
            m1_to_medoid = float(P[0, med])
            m1_centrality = float(P[0].mean())            # model-1 mean dist to others
            centrality = P.mean(1)
            m1_pct = float((centrality < centrality[0]).mean())
        else:
            spread = 0.0; med = 0; m1_to_medoid = 0.0; m1_centrality = 0.0; m1_pct = 0.0

        emit_all = I.kabsch_rmsd_batch(M, emit)
        # pool best against each model
        Wp = u["W"][p]
        pb = np.array([I.kabsch_rmsd_batch(Wp, mm).min() for mm in M])

        out.append(dict(
            pdb=pdb, n=n, fail18=pdb in I.FAIL18, n_models=m,
            model1_check=d_check,
            ens_spread=spread, m1_to_medoid=m1_to_medoid,
            m1_centrality=m1_centrality, m1_centrality_pct=m1_pct,
            emit_model1=float(emit_all[0]), emit_best_model=float(emit_all.min()),
            emit_mean_model=float(emit_all.mean()),
            poolbest_model1=float(pb[0]), poolbest_best_model=float(pb.min()),
        ))
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/126 free={I.free_gb():.1f}", flush=True)

    I.write("fail_ensemble", out)
    ok = [r for r in out if "error" not in r]
    f = [r for r in ok if r["fail18"]]; o = [r for r in ok if not r["fail18"]]
    ctrl = json.load(open(os.path.join(ROOT, "s12", "results", "fail_contrast.json")))["controls"]
    cset = {c["ctrl"] for c in ctrl}
    m18 = [r for r in ok if r["pdb"] in cset]
    print(f"\nparsed {len(ok)}/126;  model1 check max = "
          f"{max(r['model1_check'] for r in ok):.4f} A (must be ~0)")
    print(f"{'quantity':24s} {'FAIL18':>9s} {'MATCH18':>9s} {'other108':>9s}")
    for k in ["n_models", "ens_spread", "m1_to_medoid", "m1_centrality_pct",
              "emit_model1", "emit_best_model", "emit_mean_model",
              "poolbest_model1", "poolbest_best_model"]:
        print(f"{k:24s} {np.mean([r[k] for r in f]):9.3f} "
              f"{np.mean([r[k] for r in m18]):9.3f} {np.mean([r[k] for r in o]):9.3f}")
    print("\nper-target FAIL18 (declared metric first):")
    print(f"{'pdb':6s} {'mdl':>4s} {'spread':>7s} {'m1-med':>7s} {'m1pct':>6s} "
          f"{'emit@m1':>8s} {'emit@best':>9s} {'gain':>6s} {'pb@m1':>6s} {'pb@best':>7s}")
    for r in sorted(f, key=lambda r: r["pdb"]):
        print(f"{r['pdb']:6s} {r['n_models']:4d} {r['ens_spread']:7.3f} {r['m1_to_medoid']:7.3f} "
              f"{r['m1_centrality_pct']:6.2f} {r['emit_model1']:8.3f} {r['emit_best_model']:9.3f} "
              f"{r['emit_model1']-r['emit_best_model']:6.3f} {r['poolbest_model1']:6.3f} "
              f"{r['poolbest_best_model']:7.3f}")
    a = np.array([r["emit_model1"] - r["emit_best_model"] for r in f])
    b = np.array([r["emit_model1"] - r["emit_best_model"] for r in m18])
    print(f"\nensemble-aware gain: FAIL18 {a.mean():.3f}  MATCH18 {b.mean():.3f}")
    print(json.dumps(I.paired(a, b), indent=1))


if __name__ == "__main__":
    main()
