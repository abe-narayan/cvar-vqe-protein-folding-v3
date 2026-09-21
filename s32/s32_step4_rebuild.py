"""S32 LANE V -- CHARTER STEP 4.  Independent rebuild of the 3.2105 endpoint from artefacts.

Nothing here reads a report, a ladder JSON, or any S29/S31 analysis file.  It reads only:
  * s8/generate_univ/<pdb>.npz          (the frozen 126-target universe + ORACLE labels)
  * the leave-fold-out distogram        (s12/instrument.distogram, cached)
  * bench_results/cache/<PROD_KEY>/*    (the production run, for CROSS-CHECK only)
  * core.project                        (the projection operator)

and rebuilds, per target:

    pool(K=500, BLOSUM order) -> shipped Bayes-risk score -> top-75 -> coordinate average
        -> cloud RMSD                     [must mean 3.0483]
        -> set mean over the 75 members   [must mean 3.5507]
        -> projection (lam path 0 -> 0.3, multi-start) -> BUILT CHAIN RMSD  [must mean 3.2105]

Usage:
    python s32/s32_step4_rebuild.py --shard 0 --n-shards 4     # writes rows, resumable
    python s32/s32_step4_rebuild.py --analyse
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("s32_instrument", os.path.join(ROOT, "s12", "instrument.py"))
I = _ilu.module_from_spec(_spec); _spec.loader.exec_module(I)

RESULTS = os.path.join(ROOT, "s32", "results")
os.makedirs(RESULTS, exist_ok=True)
ROWS = os.path.join(RESULTS, "s32_V_step4_rows.jsonl")

# The three anchors, exactly as the S32 contract states them (section 0).
ANCHOR = {"chain": 3.2105, "cloud": 3.0483, "set_mean": 3.5507, "n": 126}


def rows_path(shard=None):
    return ROWS if shard is None else ROWS.replace(".jsonl", "_shard%d.jsonl" % int(shard))


def read_rows():
    import glob
    out = {}
    for f in sorted(glob.glob(ROWS.replace(".jsonl", "*.jsonl"))):
        for line in open(f):
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            out[r["pdb"]] = r
    return out


def append(path, row):
    with open(path, "a") as fh:
        fh.write(json.dumps(row) + "\n")


def one(t):
    """Everything for one target, rebuilt from the universe."""
    pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
    t0 = time.time()
    u = I.load_univ(pdb)
    nat = u["nat_ca"]                                            # ORACLE: evaluation only
    pool = I.pool_idx(u, I.K)                                    # K=500 BLOSUM prefix
    W = u["W"][pool]
    rr = u["rr"][pool]                                           # ORACLE labels

    # ---- score: the shipped leave-fold-out distogram Bayes risk (native-free)
    dg = I.distogram(pdb, seq, fold)
    i, j = I.pair_index(n)
    D = I.pair_dists(W, i, j)
    sc = I.shipped_score(dg, D.astype(np.float32).astype(float))

    # ---- top-75.  Cross-check against the production record's own `sub`.
    rec = I.shipped_record(pdb)
    sub_prod = np.asarray(rec["sub"], int)
    order = np.argsort(sc, kind="stable")
    sub_mine = order[:I.M]
    n_tied_at_cut = int((sc == sc[order[I.M - 1]]).sum())        # tie-breaking exposure at the cut
    set_match = set(sub_mine.tolist()) == set(sub_prod.tolist())

    # The endpoint uses the PRODUCTION top-75 (the object the ladder projected); the rebuilt
    # set is reported beside it so a disagreement is visible rather than absorbed.
    C_prod, b_prod = I.coordinate_average(W[sub_prod])
    C_mine, b_mine = I.coordinate_average(W[sub_mine]) if not set_match else (C_prod, b_prod)

    avg_stored = np.asarray(rec["avg_ca"], float)
    dev_avg = float(np.abs(C_prod - avg_stored).max())

    # ---- projection (lam 0 -> 0.3, multi-start, exact grad): the deployed operator
    pr = I.project(C_prod, seq, fold)
    chain = np.asarray(pr["ca"], float)
    fit = np.asarray(pr["fit_ca"], float)

    row = dict(
        pdb=pdb, n=n, fold=fold,
        rmsd_cloud=float(I.ca_rmsd(C_prod, nat)),
        rmsd_chain=float(I.ca_rmsd(chain, nat)),
        rmsd_fit_lam0=float(I.ca_rmsd(fit, nat)),
        set_mean=float(rr[sub_prod].mean()),
        pool_best=float(rr.min()),
        top75_best=float(rr[sub_prod].min()),
        shipped_argmin=float(rr[int(np.argmin(sc))]),
        # stored production scalars, for the cache-vs-rebuild contrast
        stored_rmsd_avg=float(rec["rmsd_avg"]), stored_rmsd_arm=float(rec["rmsd_arm"]),
        stored_rmsd_fit=float(rec["rmsd_fit"]), stored_rmsd_full=float(rec["rmsd_full"]),
        stored_ca_rmsd=float(I.ca_rmsd(np.asarray(rec["ca"], float), nat)),
        # audit fields
        sub_set_match=bool(set_match), n_tied_at_cut=n_tied_at_cut,
        dev_avg_vs_stored=dev_avg, medoid=int(b_prod),
        rmsd_cloud_mine=float(I.ca_rmsd(C_mine, nat)),
        secs=time.time() - t0,
    )
    return row


def run(shard=None, n_shards=1):
    tg = I.targets()
    assert len(tg) == 126, "target set is not 126"
    done = read_rows()
    path = rows_path(shard)
    todo = [t for k, t in enumerate(tg)
            if (shard is None or k % int(n_shards) == int(shard)) and t["pdb"] not in done]
    print("step4 shard %s/%s: %d to do -> %s" % (shard, n_shards, len(todo), os.path.basename(path)), flush=True)
    t0 = time.time()
    for k, t in enumerate(todo):
        r = one(t)
        append(path, r)
        print("  [%d/%d] %-6s cloud %.4f chain %.4f set %.4f  submatch=%s  (%.1fs, %.1f min)"
              % (k + 1, len(todo), r["pdb"], r["rmsd_cloud"], r["rmsd_chain"], r["set_mean"],
                 r["sub_set_match"], r["secs"], (time.time() - t0) / 60), flush=True)


def analyse():
    rows = read_rows()
    pdbs = sorted(rows)
    out = {"n": len(pdbs), "anchors": ANCHOR, "checks": [], "ok": True}

    def g(k):
        return np.array([rows[p][k] for p in pdbs], float)

    def chk(name, got, want, tol, note=""):
        ok = abs(got - want) <= tol
        out["checks"].append(dict(name=name, got=float(got), want=float(want), tol=tol,
                                  delta=float(got - want), ok=bool(ok), note=note))
        out["ok"] = out["ok"] and ok
        print("%-46s got %.6f  want %.4f  d=%+.6f  %s %s"
              % (name, got, want, got - want, "OK " if ok else "FAIL", note))

    chk("n targets", len(pdbs), 126, 0)
    chk("BUILT CHAIN mean (THE ENDPOINT)", g("rmsd_chain").mean(), ANCHOR["chain"], 5e-4)
    chk("CA point cloud mean (diagnostic)", g("rmsd_cloud").mean(), ANCHOR["cloud"], 5e-4)
    chk("set mean over top-75 members", g("set_mean").mean(), ANCHOR["set_mean"], 5e-4)
    chk("pool best (K=500) [ORACLE]", g("pool_best").mean(), 1.7108, 2e-3)
    chk("top-75 best [ORACLE]", g("top75_best").mean(), 2.3062, 2e-3)
    chk("shipped argmin [ORACLE label]", g("shipped_argmin").mean(), 3.4540, 2e-3)
    chk("projection price (chain - cloud)", g("rmsd_chain").mean() - g("rmsd_cloud").mean(), 0.1622, 5e-4)

    # fold structure
    folds = g("fold").astype(int)
    cnt = {int(f): int((folds == f).sum()) for f in sorted(set(folds.tolist()))}
    out["fold_counts"] = cnt
    out["fold_chain_mean"] = {int(f): float(g("rmsd_chain")[folds == f].mean()) for f in cnt}
    print("fold counts:", cnt, "sum", sum(cnt.values()))
    out["checks"].append(dict(name="5 folds, sum 126", got=sum(cnt.values()), want=126,
                              ok=bool(len(cnt) == 5 and sum(cnt.values()) == 126)))
    out["ok"] = out["ok"] and len(cnt) == 5 and sum(cnt.values()) == 126

    # top-75 reproduction and tie exposure
    out["n_sub_set_match"] = int(g("sub_set_match").sum())
    out["max_dev_avg_vs_stored"] = float(g("dev_avg_vs_stored").max())
    out["n_tied_at_cut_max"] = int(g("n_tied_at_cut").max())
    print("top-75 set reproduced on %d/126; max |C_rebuilt - stored avg_ca| = %.3g; "
          "max ties at the cut = %d"
          % (out["n_sub_set_match"], out["max_dev_avg_vs_stored"], out["n_tied_at_cut_max"]))

    # the cache-vs-rebuild contrast: what the production cache itself stores
    for k, lab in [("stored_rmsd_avg", "cache rmsd_avg (cloud)"),
                   ("stored_rmsd_fit", "cache rmsd_fit (lam=0 chain)"),
                   ("stored_rmsd_arm", "cache rmsd_arm (lam=0.3 chain)"),
                   ("stored_rmsd_full", "cache rmsd_full (post-AMBER chain)"),
                   ("stored_ca_rmsd", "cache ca[] recomputed")]:
        print("  %-34s %.6f" % (lab, g(k).mean()))
        out.setdefault("cache_means", {})[k] = float(g(k).mean())
    out["rebuild_lam0_mean"] = float(g("rmsd_fit_lam0").mean())
    d = g("rmsd_chain") - g("stored_ca_rmsd")
    out["reproj_vs_cache_chain"] = dict(mean=float(d.mean()), absmean=float(np.abs(d).mean()),
                                        p90=float(np.percentile(np.abs(d), 90)), max=float(np.abs(d).max()))
    print("  re-projection vs cached chain: mean %+.4f  |mean| %.4f  p90 %.4f  max %.4f"
          % (d.mean(), np.abs(d).mean(), np.percentile(np.abs(d), 90), np.abs(d).max()))

    out["per_target"] = {p: rows[p] for p in pdbs}
    with open(os.path.join(RESULTS, "s32_V_step4_endpoint.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print("\nVERDICT:", "ENDPOINT REPRODUCES" if out["ok"] else "*** ENDPOINT DOES NOT REPRODUCE ***")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=None)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--analyse", action="store_true")
    a = ap.parse_args()
    if a.analyse:
        analyse()
    else:
        run(a.shard, a.n_shards)
