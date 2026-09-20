#!/usr/bin/env python
"""s29/s29_M_harness.py -- lane M's evaluation-harness audit (S29 brief deliverable 3).

Nine checks, each printing its own number and PASS / FAIL:

  1  ca_rmsd is the Kabsch CA RMSD          (5-line independent implementation, 10 targets)
  2  the native is the pinned one           (universe nat_ca vs peptide_db vs a fresh PDB parse)
  3  no benchmark file is read on the path  (import-closure grep + an in-process POISON that makes
                                             any read of results/benchmark_manifest.json raise)
  4  pinned_folds are the pinned folds      (ST.pinned_folds vs the universes vs peptide_folds.json)
     and every leave-fold-out model in the deployed path really excludes its fold
  5  the 3.2126 anchor reproduces           (from chain_rows.jsonl :: DIS, and from a FRESH
                                             re-projection on 6 targets, under the poison of check 3)
  6  ST.compare's MDE is 2.8016 x SE        (read the code, then reproduce both numbers)
     and its fold CI is a cluster bootstrap over the 5 folds
  7  the 126 targets are the 126            (count, no duplicates, lengths 9 to 16)

Run:  python s26/jobrun.py --agent S29M --tag CPU --name m_harness_audit --est-ram 1.0 -- \
          python s29/s29_M_harness.py
Artefact: s29/results/s29_M_harness_audit.json
"""
from __future__ import annotations

import builtins
import glob
import io
import json
import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = os.path.join(RESULTS, "s29_M_harness_audit.json")

#: every artefact that would open the sealed benchmark
FORBIDDEN = ("benchmark_manifest.json", "final_report.json", "monomer_manifest.json")


# ============================================================ the poison
class BenchmarkPoison:
    """Make ANY read of a sealed-benchmark artefact raise, in this process.

    Stronger than the brief's "rename it in a temp copy": a rename is defeated by a cached
    read, and copying the tree would copy 1.1 GB of pinned caches. This intercepts
    `builtins.open`, `os.path.exists` and `np.load`, so a module that already holds the path,
    or that only probes for existence, is caught too. The pinned file is never touched --
    contract rule 2 forbids moving a pinned artefact, and this is strictly stronger.
    """

    def __init__(self):
        self.hits = []

    def __enter__(self):
        self._open, self._exists, self._load = builtins.open, os.path.exists, np.load

        def guard(path, where):
            s = str(path)
            if any(f in s for f in FORBIDDEN):
                self.hits.append((where, s))
                raise RuntimeError("BENCHMARK POISON: %s opened %s" % (where, s))

        def open_(file, *a, **k):
            guard(file, "open")
            return self._open(file, *a, **k)

        def exists_(path):
            guard(path, "os.path.exists")
            return self._exists(path)

        def load_(file, *a, **k):
            guard(file, "np.load")
            return self._load(file, *a, **k)

        builtins.open, os.path.exists, np.load = open_, exists_, load_
        return self

    def __exit__(self, *exc):
        builtins.open, os.path.exists, np.load = self._open, self._exists, self._load
        return False


# ============================================================ check 1
def ind_rmsd(a, b):
    """An independent Kabsch CA-RMSD, written from the definition, 5 lines."""
    a = np.asarray(a, float) - np.asarray(a, float).mean(0)
    b = np.asarray(b, float) - np.asarray(b, float).mean(0)
    U, S, Vt = np.linalg.svd(a.T @ b)
    R = U @ np.diag([1.0, 1.0, np.sign(np.linalg.det(U @ Vt))]) @ Vt
    return float(np.sqrt(((a @ R - b) ** 2).sum() / len(a)))


def check1(I):
    worst = 0.0
    rows = []
    for t in I.targets()[:10]:
        u = I.load_univ(t["pdb"])
        p = I.pool_idx(u)
        W, nat = u["W"][p], u["nat_ca"]
        for k in (0, 1, 7, 100, 499):
            worst = max(worst, abs(I.ca_rmsd(W[k], nat) - ind_rmsd(W[k], nat)))
        C, _ = I.coordinate_average(W[:75])
        worst = max(worst, abs(I.ca_rmsd(C, nat) - ind_rmsd(C, nat)))
        rng = np.random.default_rng(0)
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(Q) < 0:
            Q[:, 0] *= -1.0
        worst = max(worst, abs(I.ca_rmsd(W[0] @ Q.T + 17.0, nat) - I.ca_rmsd(W[0], nat)))
        mir = W[0] * np.array([1.0, 1.0, -1.0])
        rows.append(dict(pdb=t["pdb"], ca_rmsd=I.ca_rmsd(W[0], nat), independent=ind_rmsd(W[0], nat),
                         mirrored=I.ca_rmsd(mir, nat)))
    ok = worst < 1e-9
    print("[1] worst |I.ca_rmsd - independent Kabsch| over 10 targets x 6 pairs "
          "+ rigid-motion invariance: %.3e   %s" % (worst, "PASS" if ok else "FAIL"))
    return dict(name="ca_rmsd is Kabsch", pass_=bool(ok), worst_abs_diff=float(worst), rows=rows)


# ============================================================ check 2
def check2(I):
    from core import data, geometry as geo
    byp = {p.pdb: p for p in data.load()}
    worst_db, miss = 0.0, []
    for t in I.targets():
        u = I.load_univ(t["pdb"])
        p = byp.get(t["pdb"])
        if p is None or p.seq != t["seq"]:
            miss.append(t["pdb"])
            continue
        worst_db = max(worst_db, float(np.abs(np.asarray(u["nat_ca"], float) - np.asarray(p.ca, float)).max()))
    worst_parse, nparse = 0.0, 0
    for t in I.targets()[:10]:
        for d in ("pdbs", "pdbs_ext"):
            f = os.path.join(ROOT, d, t["pdb"] + ".pdb")
            if os.path.exists(f):
                seq, coords, _, _ = geo.native_coords_from_pdb(f)
                if seq == t["seq"]:
                    u = I.load_univ(t["pdb"])
                    worst_parse = max(worst_parse, float(np.abs(np.asarray(u["nat_ca"], float)
                                                               - np.asarray(coords["CA"], float)).max()))
                    nparse += 1
                break
    cat = os.path.join(ROOT, "catrace_prior.npz")
    catkeys = list(np.load(cat, allow_pickle=True).files) if os.path.exists(cat) else []
    ok = (not miss) and worst_db < 1e-9 and worst_parse < 1e-9 and nparse == 10
    print("[2] universe nat_ca vs peptide_db over 126: %.3e A; vs a fresh PDB parse on %d/10: %.3e A; "
          "catrace_prior.npz keys %s   %s" % (worst_db, nparse, worst_parse, catkeys,
                                              "PASS" if ok else "FAIL"))
    return dict(name="the native is the pinned one", pass_=bool(ok), worst_vs_peptide_db=float(worst_db),
                worst_vs_fresh_parse=float(worst_parse), n_parsed=int(nparse),
                missing=miss, catrace_prior_keys=catkeys)


# ============================================================ check 3 + check 5
def check3_5(I):
    """The poison test AND the anchor reproduction, in one pass (the poison must be live
    while the deployed path runs, or it proves nothing)."""
    from s24 import d_harness as H
    from s27 import run_pool as RP

    # (a) the static grep over the import closure of s27/run_vqe_chain.py
    closure = ["s27/run_vqe_chain.py", "s27/run_pool.py", "s27/ham_lib.py", "s25/phys_lib.py",
               "s24/d_harness.py", "s24/stats_lib.py", "s22/qcand_lib.py", "s16/energy_lib.py",
               "s15/seed.py", "s12/instrument.py", "core/quantum.py", "core/project.py",
               "core/predict.py", "core/geometry.py", "core/cache.py", "core/data.py"]
    pat = re.compile(r"benchmark_manifest|final_report\.json|monomer_manifest|db\.benchmark\(|\.benchmark\(\)")
    grep = {}
    for f in closure:
        p = os.path.join(ROOT, f)
        if not os.path.exists(p):
            continue
        hits = [(k + 1, ln.strip()[:120]) for k, ln in enumerate(open(p, encoding="utf-8").read().split("\n"))
                if pat.search(ln)]
        if hits:
            grep[f] = hits

    # (b) the anchor, from the stored rows
    rows = [json.loads(l) for l in open(os.path.join(ROOT, "s27", "results", "chain_rows.jsonl"),
                                        encoding="utf-8")]
    dis = [r for r in rows if r["config"] == "DIS"]
    stored_chain = float(np.mean([r["rmsd_chain"] for r in dis]))
    stored_cloud = float(np.mean([r["rmsd_cloud"] for r in dis]))

    # (c) a FRESH re-projection on 6 targets, with the poison live
    six = [t["pdb"] for t in I.targets()[:6]]
    fresh, worst = [], 0.0
    t0 = time.time()
    with BenchmarkPoison() as poison:
        for pdb in six:
            cand, ch, _ = RP.channels_for(pdb)
            E = RP.zr(ch["DIS"])
            key = RP.rng_for(pdb, "tiekey").random(cand.k)
            top = np.lexsort((key, E))[:RP.M]
            C, _ = H.readout_uniform(cand, top)
            ca = H.readout_projected(cand, C)
            got_cloud = float(I.ca_rmsd(C, cand.nat_ca))
            got_chain = float(I.ca_rmsd(ca, cand.nat_ca))
            ref = [r for r in dis if r["pdb"] == pdb][0]
            d_cloud = abs(got_cloud - ref["rmsd_cloud"])
            d_chain = abs(got_chain - ref["rmsd_chain"])
            worst = max(worst, d_cloud, d_chain)
            fresh.append(dict(pdb=pdb, cloud=got_cloud, cloud_row=ref["rmsd_cloud"],
                              chain=got_chain, chain_row=ref["rmsd_chain"],
                              d_cloud=d_cloud, d_chain=d_chain))
        hits = list(poison.hits)
    secs = time.time() - t0

    ok3 = (not hits) and all(f == "core/data.py" or f == "s7/debias.py" for f in grep)
    ok5 = abs(stored_chain - 3.2126) < 5e-4 and abs(stored_cloud - 3.0483) < 5e-4 and worst < 1e-6
    print("[3] poison live through 6 fresh target rebuilds: %d forbidden reads; static grep hits in %s   %s"
          % (len(hits), sorted(grep) or "nothing on the path", "PASS" if ok3 else "FAIL"))
    print("[5] chain_rows DIS: chain %.4f (anchor 3.2126), cloud %.4f (anchor 3.0483), n=%d; "
          "fresh re-projection on 6 targets, worst |diff| %.3e A (%.0f s)   %s"
          % (stored_chain, stored_cloud, len(dis), worst, secs, "PASS" if ok5 else "FAIL"))
    c3 = dict(name="no benchmark file on the tuning path", pass_=bool(ok3),
              poison_hits=hits, static_grep=grep)
    c5 = dict(name="the 3.2126 anchor reproduces", pass_=bool(ok5), n_rows=len(dis),
              stored_chain=stored_chain, stored_cloud=stored_cloud,
              worst_fresh_diff=float(worst), secs=float(secs), rows=fresh)
    return c3, c5


# ============================================================ check 4
def check4(I):
    from s24 import stats_lib as ST
    from core import data, predict as dgm
    tg = I.targets()
    pf = ST.pinned_folds()
    univ = np.array([t["fold"] for t in tg], int)
    folds_json = json.load(open(os.path.join(ROOT, "peptide_folds.json"), encoding="utf-8"))
    jf = np.array([folds_json["assign"][t["seq"]] for t in tg], int)
    same1 = bool(np.array_equal(pf, univ))
    same2 = bool(np.array_equal(univ, jf))
    sub = ST.pinned_folds([t["pdb"] for t in tg[:10]])
    same3 = bool(np.array_equal(sub, univ[:10]))

    # every leave-fold-out model in the deployed path really excludes its fold:
    # (a) the checkpoint served to a target is fold_model(target.fold)
    # (b) that fold's training set excludes every peptide of that fold
    served, excl_ok, leak = [], True, []
    for t in tg[:5] + tg[-5:]:
        f = int(t["fold"])
        path = dgm._model_path(f, True, True, 0)
        entries = [p for p in data.load() if data.folds(5)[p.seq] != f]
        in_train = any(p.seq == t["seq"] for p in entries)
        frags = data.fold_fragments(f, 5)
        frag_leak = [p.pdb for p in frags if t["seq"] in p.seq or p.seq in t["seq"]]
        served.append(dict(pdb=t["pdb"], fold=f, model=os.path.basename(path),
                           model_exists=os.path.exists(path), own_seq_in_train=bool(in_train),
                           verbatim_in_fragments=frag_leak[:3]))
        excl_ok = excl_ok and (not in_train) and os.path.exists(path)
        if frag_leak:
            leak.append((t["pdb"], len(frag_leak)))
    ok = same1 and same2 and same3 and excl_ok
    print("[4] pinned_folds == universes: %s; universes == peptide_folds.json: %s; by-pdb subset: %s; "
          "10 targets served by fold_model(own fold) with own sequence absent from training: %s; "
          "verbatim fragment self-copies: %s   %s"
          % (same1, same2, same3, excl_ok, leak or "none in the 10", "PASS" if ok else "FAIL"))
    return dict(name="pinned folds, and leave-fold-out really excludes", pass_=bool(ok),
                pinned_eq_universe=same1, universe_eq_json=same2, by_pdb_ok=same3,
                exclusion_ok=bool(excl_ok), served=served, fragment_self_copies=leak)


# ============================================================ check 6
def check6():
    from s24 import stats_lib as ST
    src = open(os.path.join(ROOT, "s24", "stats_lib.py"), encoding="utf-8").read()
    has_k = "MDE_K = 2.8016" in src
    has_cluster = "rf.choice(F, len(F), replace=True)" in src
    rng = np.random.default_rng(7)
    a = rng.normal(3.2, 0.9, 126)
    b = a + rng.normal(0.05, 0.3, 126)
    folds = np.array([i % 5 for i in range(126)])
    o = ST.compare(a, b, folds=folds, label="audit")
    d = a - b
    se = float(d.std(ddof=1) / math.sqrt(126))
    ok = (has_k and has_cluster and abs(o["se"] - se) < 1e-12
          and abs(o["mde"] - 2.8016 * se) < 1e-12 and o["n_folds"] == 5
          and o["ci95_fold"] is not None)
    print("[6] MDE_K literal present: %s; fold CI resamples FOLDS: %s; reproduced SE %.6f == %.6f, "
          "MDE %.6f == 2.8016*SE %.6f; n_folds %d   %s"
          % (has_k, has_cluster, o["se"], se, o["mde"], 2.8016 * se, o["n_folds"],
             "PASS" if ok else "FAIL"))
    return dict(name="MDE = 2.8016 x SE, fold CI is a cluster bootstrap", pass_=bool(ok),
                mde_k_literal=has_k, cluster_bootstrap_over_folds=has_cluster,
                se=o["se"], se_recomputed=se, mde=o["mde"], mde_recomputed=2.8016 * se,
                ci95_fold=o["ci95_fold"], ci95_iid=o["ci95_iid"], n_folds=o["n_folds"])


# ============================================================ check 7
def check7(I):
    tg = I.targets()
    n = len(tg)
    pdbs = [t["pdb"] for t in tg]
    seqs = [t["seq"] for t in tg]
    lens = np.array([t["n"] for t in tg], int)
    files = sorted(glob.glob(os.path.join(ROOT, "s8", "generate_univ", "*.npz")))
    ok = (n == 126 and len(set(pdbs)) == 126 and len(set(seqs)) == 126
          and lens.min() >= 9 and lens.max() <= 16 and len(files) == 126
          and pdbs == sorted(pdbs))
    print("[7] %d targets, %d distinct pdb ids, %d distinct sequences, lengths %d..%d, "
          "%d universe files, pdb-sorted: %s   %s"
          % (n, len(set(pdbs)), len(set(seqs)), lens.min(), lens.max(), len(files),
             pdbs == sorted(pdbs), "PASS" if ok else "FAIL"))
    return dict(name="the 126 are the 126", pass_=bool(ok), n=int(n), n_distinct_pdb=len(set(pdbs)),
                n_distinct_seq=len(set(seqs)), len_min=int(lens.min()), len_max=int(lens.max()),
                n_universe_files=len(files), length_histogram={int(k): int(v) for k, v in
                                                               zip(*np.unique(lens, return_counts=True))})


# ============================================================ check 8
def check8(I, n_targets=126):
    """A COLD distogram recomputation on every target, with the poison live: exercises the MLP
    path the s12 cache hides, and asks whether the cache is stale IN A WAY THAT MATTERS.

    Two separate questions, reported separately because they have different answers:
      (a) is the cached posterior BIT-IDENTICAL to a fresh recomputation?  -- NO, see below
      (b) does the difference change the score ORDER, the top-75 SET, or the emitted RMSD? -- the
          pass criterion, because that is the only way a stale cache could move a number.
    """
    from core import pipeline as pl, predict as dgm
    from s24 import d_harness as H
    from s27 import run_pool as RP
    tg = I.targets()[:n_targets]
    pl.guard_esm([t["seq"] for t in tg])
    worst_prob = worst_risk = worst_rel = worst_score = worst_rmsd = 0.0
    order_diff = set_diff = 0
    bad = []
    t0 = time.time()
    with BenchmarkPoison() as poison:
        # determinism first: two calls in ONE process must agree bit-for-bit
        model0 = pl.fold_model(int(tg[0]["fold"]))
        a = dgm.Distogram.for_target(tg[0]["seq"], model=model0)
        b = dgm.Distogram.for_target(tg[0]["seq"], model=model0)
        in_process = float(np.abs(a.prob - b.prob).max())
        for t in tg:
            cached = I.distogram(t["pdb"], t["seq"], t["fold"])
            d = dgm.Distogram.for_target(t["seq"], model=pl.fold_model(int(t["fold"])))
            wp = float(np.abs(np.asarray(d.prob, np.float32) - np.asarray(cached["prob"], np.float32)).max())
            wr = float(np.abs(np.asarray(d._risk, np.float32) - np.asarray(cached["risk"], np.float32)).max())
            rel = float((np.abs(np.asarray(d._risk, np.float32) - np.asarray(cached["risk"], np.float32))
                         / np.maximum(np.abs(np.asarray(cached["risk"], float)), 1e-9)).max())
            worst_prob, worst_risk, worst_rel = max(worst_prob, wp), max(worst_risk, wr), max(worst_rel, rel)
            cand = H.Candidates.from_universe(t["pdb"], k=500)
            i, j = I.pair_index(cand.n)
            D = I.pair_dists(cand.W, i, j).astype(np.float32).astype(float)
            sc_c = I.shipped_score(cached, D)
            sc_f = I.shipped_score({"grid": np.asarray(d.grid, np.float32),
                                    "risk": np.asarray(d._risk, np.float32)}, D)
            worst_score = max(worst_score, float(np.abs(sc_c - sc_f).max()))
            same_order = bool(np.array_equal(np.argsort(sc_c, kind="stable"), np.argsort(sc_f, kind="stable")))
            key = RP.rng_for(t["pdb"], "tiekey").random(cand.k)
            tc = set(np.lexsort((key, RP.zr(sc_c)))[:RP.M].tolist())
            tf = set(np.lexsort((key, RP.zr(sc_f)))[:RP.M].tolist())
            same_set = tc == tf
            Cc, _ = H.readout_uniform(cand, np.array(sorted(tc), int))
            Cf, _ = H.readout_uniform(cand, np.array(sorted(tf), int))
            dr = abs(I.ca_rmsd(Cc, cand.nat_ca) - I.ca_rmsd(Cf, cand.nat_ca))
            worst_rmsd = max(worst_rmsd, float(dr))
            order_diff += (not same_order)
            set_diff += (not same_set)
            if (not same_set) or dr > 0:
                bad.append(dict(pdb=t["pdb"], same_order=same_order, same_set=same_set, d_rmsd=float(dr)))
        hits = list(poison.hits)
    ok = (not hits) and in_process == 0.0 and set_diff == 0 and worst_rmsd == 0.0
    print("[8] cold distogram on %d targets, poison live (%d forbidden reads): in-process determinism "
          "%.1e; vs the s12 cache max|prob| %.2e, max|risk| %.2e (rel %.2e), max|score| %.2e; "
          "score ORDER differs on %d/%d, top-75 SET differs on %d/%d, max |RMSD diff| %.2e A (%.0f s)   %s"
          % (len(tg), len(hits), in_process, worst_prob, worst_risk, worst_rel, worst_score,
             order_diff, len(tg), set_diff, len(tg), worst_rmsd, time.time() - t0,
             "PASS" if ok else "FAIL"))
    return dict(name="the cached posterior is not stale in any way that moves a number",
                pass_=bool(ok), n_targets=len(tg), poison_hits=hits,
                in_process_determinism=in_process, max_abs_prob_diff=worst_prob,
                max_abs_risk_diff=worst_risk, max_rel_risk_diff=worst_rel,
                max_abs_score_diff=worst_score, n_order_differs=int(order_diff),
                n_top75_set_differs=int(set_diff), max_abs_rmsd_diff=worst_rmsd, offenders=bad,
                note=("the cache is NOT bit-identical to a fresh recomputation (float32-scale, "
                      "in-process determinism is exact) -- declared, and it changes no set and no RMSD"))


# ============================================================ check 9
def check9(I):
    """EMPIRICAL leave-fold-out: model f must NOT fit its own fold's targets better than the
    other folds' targets. ORACLE DIAGNOSTIC (it reads the native distances); it tunes nothing.

    A checkpoint trained with the wrong exclusion would show a LOWER NLL on its own fold.
    """
    from core import pipeline as pl, predict as dgm
    tg = I.targets()
    seqs = [t["seq"] for t in tg]
    pl.guard_esm(seqs)
    nll = np.zeros((5, len(tg)))
    for f in range(5):
        model = pl.fold_model(f)
        for k, t in enumerate(tg):
            u = I.load_univ(t["pdb"])
            d = dgm.Distogram.for_target(t["seq"], model=model)
            nat = np.asarray(u["nat_ca"], float)
            dist = np.linalg.norm(nat[d.i] - nat[d.j], axis=1)
            b = np.digitize(dist, dgm.BIN_EDGES)
            p = np.clip(d.prob[np.arange(len(b)), b], 1e-12, 1.0)
            nll[f, k] = float(-np.log(p).mean())
    folds = np.array([t["fold"] for t in tg], int)
    own = np.array([nll[folds[k], k] for k in range(len(tg))])
    others = np.array([nll[[f for f in range(5) if f != folds[k]], k].mean() for k in range(len(tg))])
    deployed = float(own.mean())
    seen = float(others.mean())
    delta = deployed - seen                      # POSITIVE = the deployed model is worse = correct
    sd = float((own - others).std(ddof=1))
    se = sd / math.sqrt(len(tg))
    ok = delta > 0 and delta > 2.8016 * se
    print("[9] ORACLE DIAGNOSTIC, empirical leave-fold-out: mean NLL of the DEPLOYED (own-fold) "
          "model %.5f vs the four models that SAW this fold %.5f; delta %+.5f (SE %.5f, %.2fx MDE) "
          "-- positive means the deployed model is genuinely worse on its own held-out fold   %s"
          % (deployed, seen, delta, se, abs(delta) / (2.8016 * se) if se > 0 else float("nan"),
             "PASS" if ok else "FAIL"))
    return dict(name="EMPIRICAL leave-fold-out exclusion (ORACLE diagnostic)", pass_=bool(ok),
                nll_deployed_own_fold=deployed, nll_models_that_saw_the_fold=seen,
                delta=float(delta), se=float(se), over_mde=float(abs(delta) / (2.8016 * se)) if se > 0 else None,
                per_fold_own={int(f): float(own[folds == f].mean()) for f in range(5)})


def main():
    from s12 import instrument as I
    from s24 import stats_lib as ST
    t0 = time.time()
    out = {}
    out["check1_ca_rmsd"] = check1(I)
    out["check2_native"] = check2(I)
    c3, c5 = check3_5(I)
    out["check3_no_benchmark"] = c3
    out["check5_anchor"] = c5
    out["check4_folds"] = check4(I)
    out["check6_stats"] = check6()
    out["check7_targets"] = check7(I)
    out["check8_cold_distogram"] = check8(I)
    out["check9_lfo_empirical"] = check9(I)
    out["all_pass"] = bool(all(v["pass_"] for k, v in out.items() if k.startswith("check")))
    out["secs"] = time.time() - t0
    ST.save_atomic(OUT, out, module_file=__file__)
    print("\nALL CHECKS PASS: %s   (%.0f s)   -> %s" % (out["all_pass"], out["secs"], OUT))
    return 0 if out["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
