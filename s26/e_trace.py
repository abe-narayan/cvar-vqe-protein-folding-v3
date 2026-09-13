"""s26/e_trace.py -- the production dataflow for ONE target, stage by stage, with shapes.

Lane E (Examiner), Sprint 26, EXAMINATION.md section B.  Mirrors `core.pipeline.run_target`
stage for stage on the PRODUCTION config and prints every intermediate the persisted record
(`bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json`) does not store: the ESM arrays, the
distogram posterior, the K=500 pool, the Bayes-risk score, the top-128 prefix, the selector's
Hamiltonian E = zrank(sorted scores), the trained state, the realised CVaR tail, the
coordinate average and the projection.  Every stage output is compared with the stored
record where the record carries it.  AMBER is NOT run unless `--amber` is passed (tag the
job AMBER then); its values are read from the record.

Natives are opened only in the labelled REPORTING block at the end, after every structure
is final, exactly as `core.pipeline.label` does, and through `s12.instrument.ca_rmsd`,
which is a different RMSD implementation from the one that wrote the record.

    python s26/jobrun.py --agent E --tag CPU --name e_trace_<pdb> --est-ram 1.5 -- \
        python s26/e_trace.py <pdb>

Writes `s26/results/e_trace_<pdb>.json` (provenance-stamped).
"""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np                                                          # noqa: E402

PROD_KEY = "1fc9f2dcf489e2fb"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def bonds(ca):
    ca = np.asarray(ca, float)
    return np.linalg.norm(np.diff(ca, axis=0), axis=1)


def maxdiff(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.shape != b.shape:
        return f"SHAPE {a.shape} vs {b.shape}"
    return float(np.max(np.abs(a - b)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("--amber", action="store_true", help="also run stage 4 (tag AMBER)")
    ap.add_argument("--no-fresh", action="store_true",
                    help="skip the end-to-end run_target + label reproduction")
    a = ap.parse_args(argv)
    pdb = a.pdb.upper()
    t_start = time.perf_counter()

    import core
    from core import pipeline as pl
    from core import predict as dgm_core
    from core import data as cdata
    from s12 import instrument as I
    from s24 import stats_lib as ST

    out = {"pdb": pdb, "backends": core.backend_report(), "stages": {}}
    log = []

    def say(s=""):
        print(s, flush=True)
        log.append(s)

    cfg = pl.PROD
    cfg_q = dataclasses.replace(pl.PROD, quantum=True)
    say(f"PROD config key {cfg.key()}   (quantum=True key {cfg_q.key()})")
    say(f"backends: {out['backends']}")

    # ---------------------------------------------------------------- target
    db = core.backend("data")
    dgm = core.backend("predict")
    qm = core.backend("quantum")
    targets = {p.pdb: p for p in db.load()}
    target = targets[pdb]
    folds = db.folds(cfg.n_folds)
    fold = int(folds[target.seq])
    rec_path = os.path.join(ROOT, "bench_results", "cache", PROD_KEY, f"{pdb}.json")
    with open(rec_path) as fh:
        rec = json.load(fh)
    say(f"\n== TARGET {pdb}  seq {target.seq}  n={target.n}  fold={fold}  "
        f"(record fold {rec['fold']}, record n {rec['n']})")
    out["target"] = {"seq": target.seq, "n": int(target.n), "fold": fold,
                     "record_fold": rec["fold"], "record_path": os.path.relpath(rec_path, ROOT),
                     "record_cfg_key": rec["cfg_key"]}
    assert fold == rec["fold"], "fold disagrees with the persisted record"

    # ---------------------------------------------------------------- ESM
    pl.guard_esm([target.seq])
    hot = cdata._hot()
    in_hot = target.seq in hot
    emb_raw, con = cdata.esm_raw(target.seq)
    emb32 = cdata.esm_embed(target.seq)
    X, pi, pj = dgm.features(target.seq, True)
    say(f"\n== ESM-2  cache file {os.path.relpath(cdata.ESM_SMALL, ROOT)} "
        f"(seq in hot cache: {in_hot}; esm_cache.npz never opened)")
    say(f"   raw per-residue embedding {emb_raw.shape} {emb_raw.dtype}; contact map {con.shape}; "
        f"PCA-32 embedding {emb32.shape} (esm_pca.npz)")
    say(f"   distogram feature matrix X {X.shape} {X.dtype}  over {len(pi)} pairs (min_sep 2)")
    out["stages"]["esm"] = {"hot_cache": os.path.relpath(cdata.ESM_SMALL, ROOT),
                            "seq_in_hot_cache": bool(in_hot),
                            "raw_shape": list(emb_raw.shape), "contacts_shape": list(con.shape),
                            "pca_shape": list(emb32.shape), "features_shape": list(X.shape)}

    # ---------------------------------------------------------------- distogram
    clk = pl.Clock()
    model = pl.fold_model(fold, cfg.n_folds)
    mpath = dgm._model_path(fold, True, True, 0)
    d = dgm.Distogram.for_target(target.seq, model=model)
    say(f"\n== DISTOGRAM  fold model {os.path.relpath(mpath, ROOT)}  sha256 {sha256(mpath)[:16]}  "
        f"({len(model)} model(s) in the ensemble)")
    say(f"   NBINS {dgm.NBINS}; CENTRES {np.round(dgm.CENTRES, 3).tolist()}")
    say(f"   prob {d.prob.shape}; expected {d.expected.shape} [{d.expected.min():.3f}, "
        f"{d.expected.max():.3f}] A; sd [{d.sd.min():.3f}, {d.sd.max():.3f}]; "
        f"weights w [{d.w.min():.3f}, {d.w.max():.3f}] (mean 1)")
    say(f"   risk table _risk {d._risk.shape} {d._risk.dtype} on grid 2.0..40.0 by 0.05 "
        f"({len(d.grid)} points)")
    out["stages"]["distogram"] = {
        "model_file": os.path.relpath(mpath, ROOT), "model_sha256": sha256(mpath),
        "n_models": len(model), "nbins": int(dgm.NBINS),
        "centres": [float(x) for x in dgm.CENTRES], "prob_shape": list(d.prob.shape),
        "expected_min": float(d.expected.min()), "expected_max": float(d.expected.max()),
        "sd_min": float(d.sd.min()), "sd_max": float(d.sd.max()),
        "risk_shape": list(d._risk.shape)}

    # ---------------------------------------------------------------- retrieve
    pool = pl.retrieve(target, fold, cfg, clk)
    say(f"\n== RETRIEVE  universe n_windows={pool['n_windows']}  K={cfg.k}  "
        f"pool W {pool['W'].shape} (float32 round-tripped), W64 {pool['W64'].shape}; "
        f"sim [{pool['sim'].min():.0f}, {pool['sim'].max():.0f}]")
    say(f"   record n_windows {rec['n_windows']}  ->  match {rec['n_windows'] == pool['n_windows']}")
    src = pool["src"]
    n_pep = int(sum(1 for s in src if "_" not in str(s)))
    say(f"   pool provenance: {n_pep} windows from database peptides, {len(src) - n_pep} "
        f"from protein fragments; {len(set(map(str, src)))} distinct parents")
    out["stages"]["retrieve"] = {"n_windows": int(pool["n_windows"]), "k": int(cfg.k),
                                 "W_shape": list(pool["W"].shape),
                                 "sim_min": float(pool["sim"].min()),
                                 "sim_max": float(pool["sim"].max()),
                                 "n_from_peptides": n_pep,
                                 "n_distinct_parents": len(set(map(str, src))),
                                 "record_n_windows_match": bool(rec["n_windows"] == pool["n_windows"])}

    # ---------------------------------------------------------------- score
    pool = pl.score(pool, target.seq, fold, target.n, cfg, clk)
    sc = np.asarray(pool["sc"], float)
    say(f"\n== SCORE  sc {sc.shape} [{sc.min():.5f}, {sc.max():.5f}]  n_pairs {pool['n_pairs']}  "
        f"argmin {int(np.argmin(sc))}  ties at min {(sc == sc.min()).sum()}")
    out["stages"]["score"] = {"sc_shape": list(sc.shape), "sc_min": float(sc.min()),
                              "sc_max": float(sc.max()), "n_pairs": int(pool["n_pairs"]),
                              "argmin": int(np.argmin(sc)),
                              "n_distinct_scores": int(len(np.unique(sc)))}

    # ---------------------------------------------------------------- filter (PROD and quantum)
    sub, Wsub, Ps, top, Pt = pl.filter_pool(pool, cfg, clk)
    sub_q, _, _, top_q, Pt_q = pl.filter_pool(pool, cfg_q, clk)
    rec_sub = np.asarray(rec["sub"], int)
    say(f"\n== FILTER  PROD: sub (top-{cfg.m}) {sub.shape}, top {top.shape}, Pt {Pt.shape}; "
        f"quantum cfg: top {top_q.shape}, Pt {Pt_q.shape}")
    say(f"   sub == record['sub'] exactly: {np.array_equal(sub, rec_sub)}   "
        f"sub_q == sub: {np.array_equal(sub_q, sub)}   top_q[:75] == sub: "
        f"{np.array_equal(top_q[:cfg.m], sub)}")
    say(f"   sub[:10] {sub[:10].tolist()}   n_distinct {pool['n_distinct']} of n_top {pool['n_top']}")
    out["stages"]["filter"] = {"m": int(cfg.m), "sub_shape": list(sub.shape),
                               "top_shape_prod": list(top.shape), "top_shape_quantum": list(top_q.shape),
                               "Pt_shape_prod": list(Pt.shape), "Pt_shape_quantum": list(Pt_q.shape),
                               "sub_matches_record": bool(np.array_equal(sub, rec_sub)),
                               "sub_first10": sub[:10].tolist(),
                               "n_distinct": int(pool["n_distinct"]), "n_top": int(pool["n_top"])}

    # ---------------------------------------------------------------- Hamiltonian E
    dim = 1 << cfg_q.vqe_qubits
    o = np.asarray(top_q[:dim], int)
    E = pl._zrank(sc[o])
    r = np.arange(1, dim + 1, dtype=float)
    z = (r - r.mean()) / r.std()
    dev = float(np.max(np.abs(E - z)))
    rng_E = float(E.max() - E.min())
    say(f"\n== HAMILTONIAN  H = diag(E), E = _zrank(sc[top[:{dim}]])  E {E.shape}")
    say(f"   E[:5]  {np.round(E[:5], 6).tolist()}")
    say(f"   E[-5:] {np.round(E[-5:], 6).tolist()}")
    say(f"   standardised ranks 1..{dim}: z[:3] {np.round(z[:3], 6).tolist()} ... z[-1] {z[-1]:.6f}")
    say(f"   max |E - z| = {dev:.3e}  =  {100 * dev / rng_E:.3f}% of E's range {rng_E:.4f}; "
        f"ties in sc[o]: {dim - len(np.unique(sc[o]))}")
    out["stages"]["hamiltonian"] = {"dim": dim, "E_first5": E[:5].tolist(), "E_last5": E[-5:].tolist(),
                                    "E_range": rng_E, "max_dev_from_rank_ladder": dev,
                                    "dev_frac_of_range": dev / rng_E,
                                    "n_ties_in_prefix": int(dim - len(np.unique(sc[o])))}

    # ---------------------------------------------------------------- quantum stage
    qs = pl.quantum_stage(pool, top_q, Pt_q, fold, cfg_q, clk)
    alpha, T = pl.VQE_LFO[fold % len(pl.VQE_LFO)]
    p = qs["p"]
    val, q, mass = qm.cvar_from_probs(E, p, alpha)
    n_tail = int((mass > 0).sum())
    order = np.argsort(E, kind="stable")
    cum = np.cumsum(p[order])
    k_prefix = int(np.searchsorted(cum, alpha, side="left")) + 1
    circ = qs.get("_circ")
    say(f"\n== CVaR-VQE  n={cfg_q.vqe_qubits} qubits, layers={cfg_q.vqe_layers}, "
        f"params={cfg_q.vqe_qubits * cfg_q.vqe_layers}, iters={cfg_q.vqe_iters}, seed={cfg_q.vqe_seed}")
    say(f"   VQE_LFO[fold {fold}] = (alpha {alpha}, T {T});  cvar {qs['cvar']:.6f}; "
        f"entropy {qs['entropy_bits']:.4f} bits of {cfg_q.vqe_qubits}; collapsed {qs['collapsed']}")
    say(f"   p {p.shape} sums to {p.sum():.12f}; max p {p.max():.5f} at prefix position "
        f"{int(np.argmax(p))}; realised tail (states with positive CVaR mass) = {n_tail}; "
        f"prefix length at alpha = {min(k_prefix, dim)}")
    say(f"   selection: sel (pool index) {qs['sel']}, sel_uniform {qs['sel_uniform']}; "
        f"sel in top-75: {int(qs['sel']) in set(sub.tolist())}")
    out["stages"]["quantum"] = {"n": cfg_q.vqe_qubits, "layers": cfg_q.vqe_layers,
                                "params": cfg_q.vqe_qubits * cfg_q.vqe_layers,
                                "iters": cfg_q.vqe_iters, "seed": cfg_q.vqe_seed,
                                "alpha": alpha, "T": T, "cvar": qs["cvar"],
                                "entropy_bits": qs["entropy_bits"], "collapsed": qs["collapsed"],
                                "p_max": float(p.max()), "realised_tail": n_tail,
                                "prefix_len_at_alpha": int(min(k_prefix, dim)),
                                "sel": qs["sel"], "sel_uniform": qs["sel_uniform"]}

    # ---------------------------------------------------------------- average
    C, b = pl.average(Wsub, Ps, clk)
    bC = bonds(C)
    say(f"\n== AVERAGE  avg_ca {C.shape}; medoid local index {b} (pool index {int(sub[b])}); "
        f"mean virtual Ca-Ca bond {bC.mean():.4f} A, min {bC.min():.4f}, max {bC.max():.4f}")
    say(f"   max |avg_ca - record avg_ca| = {maxdiff(C, rec['avg_ca'])}")
    out["stages"]["average"] = {"avg_ca_shape": list(C.shape), "medoid_local": int(b),
                                "medoid_pool_index": int(sub[b]), "bond_mean": float(bC.mean()),
                                "bond_min": float(bC.min()), "bond_max": float(bC.max()),
                                "maxdiff_vs_record": maxdiff(C, rec["avg_ca"])}

    # ---------------------------------------------------------------- projection
    ca, phi, psi, fit = pl.project(C, target.seq, fold, cfg, clk)
    bca = bonds(ca)
    say(f"\n== PROJECTION  ca {ca.shape}, phi {phi.shape}, psi {psi.shape}, fit_ca {fit.shape}; "
        f"penalty {cfg.penalty}@{cfg.lam}, multi_start {cfg.multi_start}, maxiter {cfg.maxiter}, "
        f"grad {cfg.project_grad}")
    say(f"   built-chain virtual bond mean {bca.mean():.4f} A, min {bca.min():.4f}, max {bca.max():.4f}")
    say(f"   max |ca - record ca| = {maxdiff(ca, rec['ca'])};  |phi - rec| = {maxdiff(phi, rec['phi'])};"
        f"  |psi - rec| = {maxdiff(psi, rec['psi'])};  |fit_ca - rec| = {maxdiff(fit, rec['fit_ca'])}")
    out["stages"]["projection"] = {"ca_shape": list(ca.shape), "bond_mean": float(bca.mean()),
                                   "bond_min": float(bca.min()), "bond_max": float(bca.max()),
                                   "maxdiff_ca": maxdiff(ca, rec["ca"]),
                                   "maxdiff_phi": maxdiff(phi, rec["phi"]),
                                   "maxdiff_psi": maxdiff(psi, rec["psi"]),
                                   "maxdiff_fit": maxdiff(fit, rec["fit_ca"])}

    # ---------------------------------------------------------------- AMBER
    say(f"\n== AMBER  (from the record unless --amber)  k={cfg.amber_k}, steps={cfg.amber_steps}")
    say(f"   record: amber_e0 {rec['amber_e0']}  amber_e1 {rec['amber_e1']}  amber_moved {rec['amber_moved']}"
        f"  strain_after {rec['amber_strain_after']}  err {rec['amber_err']}")
    amb = {"record_e0": rec["amber_e0"], "record_e1": rec["amber_e1"],
           "record_moved": rec["amber_moved"], "record_strain_after": rec["amber_strain_after"]}
    if rec.get("amber_ca") is not None:
        bam = bonds(rec["amber_ca"])
        say(f"   record amber_ca bond mean {bam.mean():.4f}, min {bam.min():.4f}, max {bam.max():.4f}")
        amb["record_bond_mean"] = float(bam.mean())
    if a.amber:
        r4 = pl.relax(target.seq, phi, psi, cfg, clk)
        if "ca" in r4:
            say(f"   FRESH: e0 {r4['e0']:.4f} e1 {r4['e1']:.4f} moved {r4['moved']:.6f}; "
                f"max |amber_ca - record| = {maxdiff(r4['ca'], rec['amber_ca'])}")
            amb.update({"fresh_e0": r4["e0"], "fresh_e1": r4["e1"], "fresh_moved": r4["moved"],
                        "maxdiff_amber_ca": maxdiff(r4["ca"], rec["amber_ca"])})
        else:
            say(f"   FRESH AMBER failed: {r4.get('err')}")
            amb["fresh_err"] = r4.get("err")
    out["stages"]["amber"] = amb

    # ---------------------------------------------------------------- REPORTING (natives opened here)
    nat = pl._q(np.asarray(target.ca, float), cfg)          # the reference's float32 round trip
    rr = {"rmsd_avg": I.ca_rmsd(C, nat), "rmsd_fit": I.ca_rmsd(fit, nat),
          "rmsd_arm": I.ca_rmsd(ca, nat)}
    if rec.get("amber_ca") is not None:
        rr["rmsd_full"] = I.ca_rmsd(np.asarray(rec["amber_ca"], float), nat)
    say(f"\n== REPORTING (native opened here, after every structure is final; s12.instrument.ca_rmsd)")
    for k, v in rr.items():
        say(f"   {k:10} fresh {v:.6f}   record {rec[k]:.6f}   diff {v - rec[k]:+.3e}")
    say(f"   record shipped {rec['shipped']:.6f}  pool_best {rec['pool_best']:.6f}  "
        f"top_m_best {rec['top_m_best']:.6f}  pool_mean {rec['pool_mean']:.6f}")
    out["reporting"] = {k: {"fresh": v, "record": rec[k], "diff": v - rec[k]} for k, v in rr.items()}
    out["reporting"]["record_shipped"] = rec["shipped"]
    out["reporting"]["record_pool_best"] = rec["pool_best"]
    out["reporting"]["record_top_m_best"] = rec["top_m_best"]

    # ---------------------------------------------------------------- fresh end-to-end (no AMBER)
    if not a.no_fresh:
        cfg_na = dataclasses.replace(pl.PROD, amber=False)
        clk2 = pl.Clock()
        rec2, pool2, clk2 = pl.run_target(target, fold, cfg_na, clk2)
        lab2 = pl.label(rec2, pool2, target, clk2, cfg_na)
        say(f"\n== FRESH run_target + label (amber=False, key {cfg_na.key()})")
        for k in ("shipped", "pool_best", "top_m_best", "rmsd_avg", "rmsd_fit", "rmsd_arm"):
            say(f"   {k:10} fresh {lab2[k]:.6f}   record {rec[k]:.6f}   diff {lab2[k] - rec[k]:+.3e}")
        say(f"   fresh ca vs record ca max|diff| {maxdiff(rec2['ca'], rec['ca'])}; "
            f"sub identical {np.array_equal(np.asarray(rec2['sub'], int), rec_sub)}")
        out["fresh_run"] = {k: {"fresh": lab2[k], "record": rec[k], "diff": lab2[k] - rec[k]}
                            for k in ("shipped", "pool_best", "top_m_best", "rmsd_avg",
                                      "rmsd_fit", "rmsd_arm")}
        out["fresh_run"]["maxdiff_ca"] = maxdiff(rec2["ca"], rec["ca"])
        out["fresh_run"]["sub_identical"] = bool(np.array_equal(np.asarray(rec2["sub"], int), rec_sub))
        out["fresh_run"]["timings"] = {k: round(v, 3) for k, v in clk2.t.items() if v}

    out["timings"] = {k: round(v, 3) for k, v in clk.t.items() if v}
    out["wall_s"] = round(time.perf_counter() - t_start, 2)
    out["rss"] = pl.proc_rss()
    out["log"] = log
    os.makedirs(os.path.join(ROOT, "s26", "results"), exist_ok=True)
    path = os.path.join(ROOT, "s26", "results", f"e_trace_{pdb}.json")
    ST.save_atomic(path, out, module_file=__file__)
    say(f"\nwrote {os.path.relpath(path, ROOT)}   wall {out['wall_s']} s   rss {out['rss']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
