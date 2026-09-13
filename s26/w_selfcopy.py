#!/usr/bin/env python
"""s26/w_selfcopy.py -- the 2/60 benchmark self-copy leak, bounded from the dev-set proxy.
Lane W (Wildcard), Sprint 26.  Pre-registration: s26/PREREG_selfcopy_bound.md.

THE QUESTION.  Two sealed-benchmark targets are carried verbatim by a longer database peptide
that sits in another pinned fold (S24 L4; count re-derived by lane I, L15/L18).  The leak has two
channels: (A) the carrier's self-window sits in the K = 500 retrieval pool at BLOSUM rank 0;
(B) the carrier's native distances were training labels for the target's fold model.  Four dev
targets have the same situation: 1CEK in 1A11, 2FBU in 2LMF, 2P5H in 2P5J, 6B9K in 1U6V.  This
module measures both channels on those four, a matched control population, an n = 126 upper
envelope (a fold model that trained on the target's OWN native), and turns them into a bound
on every benchmark figure, without reading any benchmark sequence, native, RMSD or the manifest.

NATIVE-FREE / GATED SPLIT.  `census`, `retrieval`, `envelope` and `posterior` read no native
quantity: every universe is loaded through `load_blind`, which overwrites `rr` and `nat_ca`
with NaN, and every stored array is asserted finite.  They store the EMISSIONS (clouds and
chains) so that `endpoint` and `floor`, which read natives, are a few seconds of RMSD
arithmetic on stored arrays and refuse to run until "PHASE 0 SIGNED OFF" is in s26/LEDGER.md.

THE TRIANGLE BOUND.  Kabsch CA-RMSD after optimal superposition is a metric on shapes, so
|RMSD(x, nat) - RMSD(y, nat)| <= RMSD(x, y).  The RMSD between the leaked and the un-leaked
emission is therefore a rigorous per-target upper bound on the leak's effect on RMSD-to-native
that never reads the native.  (Checked numerically in s26/w_selfcopy_test.py.)

COMMANDS
    census      native-free: self-windows, BLOSUM ranks, pool position, top-75 membership, ties
    retrieval   native-free: channel A emissions with/without, triangle bounds, 122 controls,
                the >= 0.6 variant (S10-4's operator on the production basis), ORACLE insertion
    envelope    native-free: the 4 leaked pinned fold models per target, emissions and bounds
    train       channel B: retrain one fold model with one chain removed (through p_ladder)
    posterior   native-free: carrier-out / control-out / reference posteriors and emissions
    endpoint    GATED: signed RMSD-to-native deltas from the stored emissions; ST.compare for C
    floor       GATED: same sequence, different deposit (the 18 verbatim dev cases)
    report      the bound arithmetic (Part D) and its assumptions

    python s26/jobrun.py --agent W --tag CPU --name w_selfcopy_census --est-ram 0.3 -- python s26/w_selfcopy.py census
    python s26/jobrun.py --agent W --tag CPU --name w_selfcopy_retrieval_probe --est-ram 0.6 -- python s26/w_selfcopy.py retrieval --probe 1CEK
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from core import data as D                           # noqa: E402
from core import predict as PR                       # noqa: E402

RES = os.path.join(HERE, "results")
MODELS = os.path.join(HERE, "models", "w_selfcopy")
P_MODELS = os.path.join(HERE, "models", "p_ladder")
os.makedirs(RES, exist_ok=True)
os.makedirs(MODELS, exist_ok=True)

K, TOPM = 500, 75
LEDGER = os.path.join(HERE, "LEDGER.md")
SIGNOFF = "PHASE 0 SIGNED OFF"
N_BENCH, N_LEAKED_BENCH = 60, 2
BENCH_CI_HALF = 0.170                                #: S9-10: +0.0103 [-0.1596, +0.1803]
IMMATERIAL = 0.1 * BENCH_CI_HALF                     #: one tenth of the benchmark CI half-width

#: (copy, carrier, copy fold, carrier fold): the four declared dev self-copies (state brief 3,
#: S24 L4, lane I L15).  The carriers are database peptides; their ids are on the record.
SELF = (("1CEK", "1A11", 2, 3), ("2FBU", "2LMF", 4, 0), ("2P5H", "2P5J", 4, 2), ("6B9K", "1U6V", 0, 2))
SELF_PDBS = tuple(s[0] for s in SELF)
#: channel-B control chains: dev-target chains (public ids) in each fold's TRAINING set, the two
#: longest that are not verbatim relatives of any dev target; fixed here before any run.
CONTROL_CHAINS = {0: ("9BAF", "8TXS"), 2: ("8TXS", "8T63"), 4: ("9BAF", "8T63")}
#: the fixed order of channel-B training jobs (PREREG section 8)
TRAIN_ORDER = (("1A11", 2), ("2LMF", 4), ("2P5J", 4), ("1U6V", 0),
               ("9BAF", 0), ("8TXS", 0), ("8TXS", 2), ("8T63", 2), ("9BAF", 4), ("8T63", 4))


# ============================================================ helpers
def signed_off():
    try:
        return SIGNOFF in open(LEDGER, encoding="utf-8").read()
    except OSError:
        return False


def require_signoff(what):
    if not signed_off():
        raise SystemExit("PHASE GATE: '%s' is not in %s; `%s` reads natives and refuses to run."
                         % (SIGNOFF, LEDGER, what))


def blind(u):
    """Poison every native quantity of a loaded universe so that no native-free command can
    read one by accident: any use would surface as NaN in the stored arrays (asserted finite)."""
    u = dict(u)
    u["rr"] = np.full(np.shape(u["rr"]), np.nan)
    u["nat_ca"] = np.full(np.shape(u["nat_ca"]), np.nan)
    return u


def load_blind(pdb):
    return blind(I.load_univ(pdb))


def targets_by_pdb():
    return {t["pdb"]: t for t in I.targets()}


def sha(s):
    return hashlib.sha256(str(s).encode()).hexdigest()[:12]


def finite(obj):
    """True if every array/number nested in obj is finite (lists, dicts, arrays)."""
    if isinstance(obj, dict):
        return all(finite(v) for v in obj.values())
    if isinstance(obj, (list, tuple)):
        return all(finite(v) for v in obj)
    if isinstance(obj, np.ndarray):
        return bool(np.isfinite(obj).all()) if obj.dtype.kind in "fc" else True
    if isinstance(obj, (float, np.floating)):
        return bool(np.isfinite(obj))
    return True


def tolist(obj):
    if isinstance(obj, dict):
        return {k: tolist(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [tolist(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def save(name, payload, rows=None, n_expected=None, complete_keys=None):
    path = os.path.join(RES, "w_selfcopy_%s.json" % name)
    ST.save_atomic(path, tolist(payload), complete_keys=complete_keys, rows=rows,
                   n_expected=n_expected, module_file=__file__)
    return path


def load_result(name):
    path = os.path.join(RES, "w_selfcopy_%s.json" % name)
    if not os.path.exists(path):
        return None
    return json.load(open(path))


# ============================================================ the deployable emission
def score_pool(W_pool, dg, n, f32=True):
    """The shipped Bayes-risk score of every pool member.  `dg` needs `grid` and `risk`
    (the s12 cache, or a `p_ladder.risk_table`).  The float32 round trip of the distances is
    the instrument's own convention (`s12.instrument.selfcheck`)."""
    i, j = I.pair_index(n)
    Dp = I.pair_dists(np.asarray(W_pool, float), i, j)
    if f32:
        Dp = Dp.astype(np.float32).astype(float)
    return np.asarray(I.shipped_score(dg, Dp), float)


def tie_set(sc):
    sc = np.asarray(sc, float)
    return np.flatnonzero(np.isclose(sc, sc.min(), atol=1e-12, rtol=1e-9))


def emit(W_pool, dg, seq, fold, project=True, f32=True):
    """DEPLOYABLE and native-free: score the pool, take the top-75 (stable argsort), average in
    the medoid frame, project onto ideal geometry (ramah 0.3, multi-start).  Returns the
    emissions and the selection bookkeeping.  Reads nothing about the native."""
    W_pool = np.asarray(W_pool, float)
    sc = score_pool(W_pool, dg, len(seq), f32=f32)
    o = np.argsort(sc, kind="stable")
    top = o[:TOPM]
    P = I.pairwise_rmsd(W_pool[top])
    C = I.superpose_batch(W_pool[top], W_pool[top][I.medoid(P)]).mean(0)
    out = {"score": sc, "top": top, "argmin_set": tie_set(sc), "cloud": C}
    if project:
        pr = I.project(C, seq, int(fold))
        out.update({"chain": pr["ca"], "fit": pr["fit_ca"], "phi": pr["phi"], "psi": pr["psi"]})
    return out


def drop_and_refill(order, drop, k=K):
    """The pool with the dropped universe indices removed and refilled from the next-best
    BLOSUM windows (S10-4's operator, the production rule).  Preserves the order."""
    order = np.asarray(order, int)
    mask = ~np.isin(order, np.asarray(drop, int))
    return order[mask][:k]


def argmin_window_rmsd(W_a, top_a, W_b, top_b, same):
    """The native-free bound on the change of the tie-averaged `sel` endpoint: zero when the two
    argmin sets are the same windows (`same`, judged by the caller on universe identity), else
    the maximum cross-pair RMSD, since |mean_a f(a) - mean_b f(b)| <= max_{a,b} |f(a) - f(b)|
    <= max_{a,b} RMSD(w_a, w_b) for f = RMSD-to-native."""
    if same:
        return 0.0
    return float(max(I.ca_rmsd(W_a[a], W_b[b]) for a in top_a for b in top_b))


def triangle(e_with, e_without, W_with, W_without, argmin_same):
    """The native-free bounds between two emissions: cloud, chain (arm), fit, and sel."""
    out = {"tri_cloud": float(I.ca_rmsd(e_with["cloud"], e_without["cloud"]))}
    if "chain" in e_with and "chain" in e_without:
        out["tri_arm"] = float(I.ca_rmsd(e_with["chain"], e_without["chain"]))
        out["tri_fit"] = float(I.ca_rmsd(e_with["fit"], e_without["fit"]))
    out["tri_sel"] = argmin_window_rmsd(W_with, e_with["argmin_set"], W_without, e_without["argmin_set"], argmin_same)
    return out


def production_gate(pdb, e, rec=None):
    """The reproduced emission against the production cache: top-75 set and the cloud."""
    rec = rec or I.shipped_record(pdb)
    sub = set(int(x) for x in rec["sub"])
    top = set(int(x) for x in e["top"])
    return {"top75_equals_sub": top == sub, "top75_overlap_with_sub": len(top & sub) / TOPM,
            "cloud_maxabs_vs_avg_ca": float(np.abs(np.asarray(e["cloud"]) - np.asarray(rec["avg_ca"], float)).max()),
            "chain_maxabs_vs_ca": (float(np.abs(np.asarray(e["chain"]) - np.asarray(rec["ca"], float)).max())
                                   if "chain" in e else None)}


# ============================================================ posteriors from other models
def pinned_posterior(seq, fold_j):
    """The posterior of the PINNED fold model `fold_j` for `seq` (as s12.instrument.distogram
    builds the shipped one).  For fold_j != the target's own fold this model trained on the
    target's native: ORACLE by construction."""
    from core import pipeline as pl
    pl.guard_esm([seq])
    model = pl.fold_model(int(fold_j))
    d = PR.Distogram.for_target(seq, model=model)
    return np.asarray(d.prob, float), np.asarray(d.i), np.asarray(d.j)


def risk_from_posterior(seq, prob, i, j):
    import p_ladder as PL
    return PL.risk_table(seq, prob, i, j)


def posterior_stats(prob_a, prob_b):
    """Native-free differences between two posteriors on the same pairs."""
    pa, pb = np.asarray(prob_a, float), np.asarray(prob_b, float)
    c = np.asarray(PR.CENTRES, float)[None, :]
    ea, eb = (pa * c).sum(1), (pb * c).sum(1)
    return {"prob_maxabs": float(np.abs(pa - pb).max()), "prob_meanabs": float(np.abs(pa - pb).mean()),
            "expected_meanabs": float(np.abs(ea - eb).mean()), "expected_maxabs": float(np.abs(ea - eb).max())}


# ============================================================ census (native-free)
def self_windows(u, seq):
    q = D.encode(seq)
    S = np.asarray(u["S"])
    return np.flatnonzero((S == q[None, :]).all(1))


def window_identities(u, seq, threshold=D.IDENTITY_THRESHOLD):
    """Shorter-normalised identity of every universe window against the target; windows are
    target-length so longer == shorter and this is S10-4's window-level measure.  Sequences
    only."""
    S = np.asarray(u["S"])
    seqs = ["".join(D.ALPHABET[int(c)] for c in row) for row in S]
    C, lens = D.composition_matrix(seqs)
    bound = np.asarray(D.max_possible_identity_many(seq, C, lens), float)
    ident = np.zeros(len(seqs))
    cand = np.flatnonzero(bound >= threshold)
    if len(cand):
        ident[cand] = D.identity_many(seq, [seqs[k] for k in cand])
    return ident


def census(verbose=True):
    tg = I.targets()
    rows = []
    for t in tg:
        pdb, seq, n = t["pdb"], t["seq"], t["n"]
        u = load_blind(pdb)
        order = np.asarray(u["order"], int); sim = np.asarray(u["sim"], float); org = np.asarray(u["org"], bool)
        rank_of = np.empty(len(order), int); rank_of[order] = np.arange(len(order))
        rec = I.shipped_record(pdb); sub = set(int(x) for x in rec["sub"])
        pool = order[:K]; pos_of = {int(w): k for k, w in enumerate(pool)}
        exact = self_windows(u, seq)
        ident = window_identities(u, seq)
        ge06 = np.flatnonzero(ident >= D.IDENTITY_THRESHOLD)
        s500 = sim[order[K - 1]]
        rank0 = int(order[0])
        row = {"pdb": pdb, "n": n, "fold": t["fold"], "self_copy": pdb in SELF_PDBS,
               "n_windows": int(len(order)), "n_exact": int(len(exact)),
               "exact_idx": exact.tolist(), "exact_blosum_rank": rank_of[exact].tolist(),
               "exact_sim": sim[exact].tolist(), "max_sim": float(sim.max()),
               "exact_from_peptide": org[exact].tolist(),
               "exact_pool_pos": [pos_of.get(int(e), -1) for e in exact],
               "exact_in_top75": [bool(pos_of.get(int(e), -1) in sub) for e in exact],
               "n_ge06_windows": int(len(ge06)), "n_ge06_in_pool": int(sum(int(g) in pos_of for g in ge06)),
               "n_ge06_in_top75": int(sum(pos_of.get(int(g), -1) in sub for g in ge06)),
               "max_identity": float(ident.max()),
               "rank0_idx": rank0, "rank0_identity": float(ident[rank0]), "rank0_in_top75": bool(0 in sub),
               "rank0_from_peptide": bool(org[rank0]),
               "boundary_sim": float(s500), "n_tied_at_boundary": int((sim == s500).sum()),
               "n_tied_inside_pool": int((sim[pool] == s500).sum())}
        rows.append(row)
        if verbose:
            print("  %s n=%2d f%d exact=%d ranks=%s pool_pos=%s in_top75=%s ge06=%d (pool %d, top75 %d) rank0_id=%.3f ties=%d/%d"
                  % (pdb, n, t["fold"], row["n_exact"], row["exact_blosum_rank"], row["exact_pool_pos"],
                     row["exact_in_top75"], row["n_ge06_windows"], row["n_ge06_in_pool"], row["n_ge06_in_top75"],
                     row["rank0_identity"], row["n_tied_inside_pool"], row["n_tied_at_boundary"]), flush=True)
        del u
    assert finite(rows)
    summ = {"n_targets": len(rows),
            "targets_with_exact_window": [r["pdb"] for r in rows if r["n_exact"]],
            "targets_with_ge06_window": [r["pdb"] for r in rows if r["n_ge06_windows"]],
            "n_targets_with_ge06_window": int(sum(1 for r in rows if r["n_ge06_windows"])),
            "exact_in_top75": {r["pdb"]: r["exact_in_top75"] for r in rows if r["n_exact"]},
            "rank0_in_top75_frac": float(np.mean([r["rank0_in_top75"] for r in rows])),
            "ties_at_boundary_median": float(np.median([r["n_tied_at_boundary"] for r in rows])),
            "ties_inside_pool_median": float(np.median([r["n_tied_inside_pool"] for r in rows])),
            "ties_inside_pool_mean": float(np.mean([r["n_tied_inside_pool"] for r in rows])),
            "targets_with_a_boundary_tie": int(sum(1 for r in rows if r["n_tied_at_boundary"] > 1))}
    p = save("census", {"label": "self-copy census, native-free", "summary": summ, "rows": rows},
             rows=rows, n_expected=len(tg), complete_keys=("pdb", "n_exact", "n_ge06_windows"))
    print(json.dumps(summ, indent=1)); print("  ->", p)
    return rows


# ============================================================ Part A: retrieval (native-free)
def same_fold_carriers(seq, fold):
    """Verbatim carriers of `seq` in the target's OWN fold (withheld by the fold discipline):
    (record, start).  Described by length only in outputs unless the carrier is a dev target."""
    folds = D.folds(5)
    out = []
    for p in D.load():
        if p.seq != seq and seq in p.seq and folds[p.seq] == fold:
            out.append((p, p.seq.find(seq)))
    return out


def retrieval_target(t, ident_threshold=D.IDENTITY_THRESHOLD, verbose=True):
    pdb, seq, n, fold = t["pdb"], t["seq"], t["n"], t["fold"]
    u = load_blind(pdb)
    W = u["W"]; order = np.asarray(u["order"], int)
    dg = I.distogram(pdb)
    rec = I.shipped_record(pdb)
    pool0 = order[:K]
    e_with = emit(W[pool0], dg, seq, fold)
    gate = production_gate(pdb, e_with, rec)
    exact = self_windows(u, seq)
    ident = window_identities(u, seq)
    ge06 = np.flatnonzero(ident >= ident_threshold)
    is_self = pdb in SELF_PDBS
    #: the dropped slot: the exact self-windows on a self-copy target, else the BLOSUM rank-0
    #: window (the matched control: same operator, an unrelated window in the slot)
    drop = exact if is_self else order[:1]
    assert (len(exact) > 0) == is_self, (pdb, len(exact), is_self)
    pool1 = drop_and_refill(order, drop)
    e_wo = emit(W[pool1], dg, seq, fold)
    refill = [int(w) for w in pool1 if w not in set(pool0.tolist())]
    argmin_same = bool(set(pool0[e_with["argmin_set"]].tolist()) == set(pool1[e_wo["argmin_set"]].tolist()))
    row = {"pdb": pdb, "n": n, "fold": fold, "self_copy": is_self, "gate": gate,
           "dropped_idx": [int(x) for x in drop], "dropped_identity": [float(ident[x]) for x in drop],
           "dropped_in_top75": [bool(int(np.flatnonzero(pool0 == x)[0]) in set(int(y) for y in e_with["top"])) for x in drop],
           "refill_idx": refill,
           "refill_in_top75": [bool(int(np.flatnonzero(pool1 == r)[0]) in set(int(y) for y in e_wo["top"])) for r in refill],
           "top75_overlap": float(len(set(pool0[e_with["top"]].tolist()) & set(pool1[e_wo["top"]].tolist())) / TOPM),
           "argmin_same": argmin_same,
           "n_tied_argmin": [int(len(e_with["argmin_set"])), int(len(e_wo["argmin_set"]))],
           **triangle(e_with, e_wo, W[pool0], W[pool1], argmin_same),
           "emissions": {"with": {k: e_with[k] for k in ("cloud", "chain", "fit")},
                         "without": {k: e_wo[k] for k in ("cloud", "chain", "fit")}},
           "argmin_windows": {"with": W[pool0][e_with["argmin_set"]], "without": W[pool1][e_wo["argmin_set"]]}}
    #: S10-4's operator on the production basis (the C27 re-derivation): drop ALL >= 0.6 windows
    if len(ge06):
        pool2 = drop_and_refill(order, ge06)
        e_ge = emit(W[pool2], dg, seq, fold)
        same_ge = bool(set(pool0[e_with["argmin_set"]].tolist()) == set(pool2[e_ge["argmin_set"]].tolist()))
        row["ge06"] = {"n_dropped": int(len(ge06)), "n_dropped_in_pool": int(np.isin(ge06, pool0).sum()),
                       "argmin_same": same_ge,
                       **triangle(e_with, e_ge, W[pool0], W[pool2], same_ge),
                       "top75_overlap": float(len(set(pool0[e_with["top"]].tolist()) & set(pool2[e_ge["top"]].tolist())) / TOPM),
                       "emissions": {k: e_ge[k] for k in ("cloud", "chain", "fit")},
                       "argmin_windows": W[pool2][e_ge["argmin_set"]]}
    #: ORACLE INSERTION: the same-fold verbatim carrier's window put into the pool (slot 500)
    ins = []
    for p, start in same_fold_carriers(seq, fold):
        win = np.asarray(p.ca[start:start + n], float)
        W_ins = np.concatenate([win[None], W[pool0][:-1]], 0)
        e_in = emit(W_ins, dg, seq, fold)
        #: W_ins index k >= 1 is pool0 position k - 1; index 0 is the inserted window
        same_in = bool(set((int(y) + 1) for y in e_with["argmin_set"] if int(y) + 1 < K) == set(int(y) for y in e_in["argmin_set"]))
        ins.append({"carrier_len": len(p.seq), "carrier_hash": sha(p.seq), "start": int(start),
                    "carrier_is_dev": p.pdb in targets_by_pdb(),
                    "inserted_in_top75": bool(0 in set(int(y) for y in e_in["top"])),
                    "inserted_score_rank": int(np.flatnonzero(np.argsort(e_in["score"], kind="stable") == 0)[0]),
                    "inserted_is_argmin": bool(0 in set(int(y) for y in e_in["argmin_set"])),
                    "argmin_same": same_in,
                    **triangle(e_with, e_in, W[pool0], W_ins, same_in),
                    "emissions": {k: e_in[k] for k in ("cloud", "chain", "fit")},
                    "argmin_windows": W_ins[e_in["argmin_set"]], "inserted_window": win})
    row["oracle_insertion"] = ins
    assert finite(row), pdb
    if verbose:
        print("  %s %s drop=%s in_top75=%s refill_in=%s | tri cloud %.3f arm %.3f fit %.3f sel %.3f | ge06 n=%d tri_arm %s | ins %d"
              % (pdb, "SELF" if is_self else "ctrl", row["dropped_idx"], row["dropped_in_top75"], row["refill_in_top75"],
                 row["tri_cloud"], row["tri_arm"], row["tri_fit"], row["tri_sel"], len(ge06),
                 ("%.3f" % row["ge06"]["tri_arm"]) if len(ge06) else "-", len(ins)), flush=True)
    del u
    return row


def retrieval(probe=None, verbose=True):
    tg = I.targets()
    if probe:
        tg = [t for t in tg if t["pdb"] == probe]
        row = retrieval_target(tg[0], verbose=verbose)
        p = save("retrieval_probe_%s" % probe, {"label": "retrieval probe", "rows": [row]})
        print("  ->", p); return [row]
    name = "retrieval"
    prev = load_result(name)
    rows = prev["rows"] if prev else []
    done = {r["pdb"] for r in rows}
    t_last = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        rows.append(retrieval_target(t, verbose=verbose))
        if time.time() - t_last > 300 or len(rows) % 10 == 0:
            save(name, {"label": "channel A, native-free", "rows": rows}, rows=rows, n_expected=len(tg))
            t_last = time.time()
    gates = [r["gate"]["top75_equals_sub"] for r in rows]
    summ = {"n": len(rows), "top75_equals_sub_count": int(sum(gates)),
            "cloud_maxabs_vs_avg_ca_max": float(max(r["gate"]["cloud_maxabs_vs_avg_ca"] for r in rows)),
            "chain_maxabs_vs_ca_max": float(max(r["gate"]["chain_maxabs_vs_ca"] for r in rows)),
            "self": {r["pdb"]: {k: r[k] for k in ("dropped_in_top75", "refill_in_top75", "top75_overlap", "argmin_same",
                                                    "tri_cloud", "tri_arm", "tri_fit", "tri_sel")}
                     for r in rows if r["self_copy"]},
            "control_tri_arm_pctiles": {q: float(np.percentile([r["tri_arm"] for r in rows if not r["self_copy"]], q))
                                        for q in (50, 90, 95, 100)},
            "control_rank0_in_top75_frac": float(np.mean([r["dropped_in_top75"][0] for r in rows if not r["self_copy"]])),
            "n_with_ge06": int(sum(1 for r in rows if "ge06" in r)),
            "n_insertions": int(sum(len(r["oracle_insertion"]) for r in rows)),
            "insertions_in_top75": int(sum(x["inserted_in_top75"] for r in rows for x in r["oracle_insertion"]))}
    p = save(name, {"label": "channel A, native-free", "summary": summ, "rows": rows}, rows=rows,
             n_expected=len(tg), complete_keys=("pdb", "tri_arm", "tri_cloud", "emissions"))
    print(json.dumps(summ, indent=1)); print("  ->", p)
    return rows


# ============================================================ Part C: envelope (native-free)
def envelope_target(t, verbose=True):
    pdb, seq, n, fold = t["pdb"], t["seq"], t["n"], t["fold"]
    u = load_blind(pdb)
    W = u["W"]; order = np.asarray(u["order"], int); pool = order[:K]; Wp = W[pool]
    dg = I.distogram(pdb)
    rec = I.shipped_record(pdb)
    e_clean = emit(Wp, dg, seq, fold)
    gate = production_gate(pdb, e_clean, rec)
    P0 = np.asarray(dg["prob"], float)
    row = {"pdb": pdb, "n": n, "fold": fold, "gate": gate,
           "clean": {k: e_clean[k] for k in ("cloud", "chain", "fit")},
           "clean_argmin_windows": Wp[e_clean["argmin_set"]], "leaked": {}}
    for j in range(5):
        if j == fold:
            continue
        prob, i, jj = pinned_posterior(seq, j)
        assert np.array_equal(i, np.asarray(dg["i"])) and np.array_equal(jj, np.asarray(dg["j"])), "pair index mismatch"
        rt = risk_from_posterior(seq, prob, i, jj)
        e_j = emit(Wp, rt, seq, fold)
        same_j = bool(set(e_clean["argmin_set"].tolist()) == set(e_j["argmin_set"].tolist()))
        row["leaked"][str(j)] = {
            **posterior_stats(prob, P0),
            "top75_overlap": float(len(set(e_clean["top"].tolist()) & set(e_j["top"].tolist())) / TOPM),
            "argmin_same": same_j,
            **triangle(e_clean, e_j, Wp, Wp, same_j),
            "emissions": {k: e_j[k] for k in ("cloud", "chain", "fit")},
            "argmin_windows": Wp[e_j["argmin_set"]]}
    assert finite(row), pdb
    if verbose:
        print("  %s f%d gate=%s | tri_arm by leaked model: %s | overlap %s"
              % (pdb, fold, gate["top75_equals_sub"],
                 " ".join("%s:%.3f" % (j, v["tri_arm"]) for j, v in row["leaked"].items()),
                 " ".join("%.2f" % v["top75_overlap"] for v in row["leaked"].values())), flush=True)
    del u
    return row


def envelope(probe=None, verbose=True):
    tg = I.targets()
    if probe:
        row = envelope_target([t for t in tg if t["pdb"] == probe][0], verbose=verbose)
        p = save("envelope_probe_%s" % probe, {"label": "envelope probe", "rows": [row]})
        print("  ->", p); return [row]
    name = "envelope"
    prev = load_result(name)
    rows = prev["rows"] if prev else []
    done = {r["pdb"] for r in rows}
    t_last = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        rows.append(envelope_target(t, verbose=verbose))
        if time.time() - t_last > 300 or len(rows) % 10 == 0:
            save(name, {"label": "Part C envelope, native-free", "rows": rows}, rows=rows, n_expected=len(tg))
            t_last = time.time()
    tri = [np.mean([v["tri_arm"] for v in r["leaked"].values()]) for r in rows]
    summ = {"n": len(rows), "top75_equals_sub_count": int(sum(r["gate"]["top75_equals_sub"] for r in rows)),
            "tri_arm_mean_over_models_pctiles": {q: float(np.percentile(tri, q)) for q in (50, 90, 100)},
            "top75_overlap_mean": float(np.mean([v["top75_overlap"] for r in rows for v in r["leaked"].values()])),
            "argmin_same_frac": float(np.mean([v["argmin_same"] for r in rows for v in r["leaked"].values()]))}
    p = save(name, {"label": "Part C envelope, native-free", "summary": summ, "rows": rows}, rows=rows,
             n_expected=len(tg), complete_keys=("pdb", "clean", "leaked"))
    print(json.dumps(summ, indent=1)); print("  ->", p)
    return rows


# ============================================================ Part B: training (channel B)
def model_path(fold, tag, seed=0):
    return os.path.join(MODELS, "pca32_fold%d_s%d_%s.pt" % (fold, seed, tag))


def reference_model_path(fold, seed=0):
    """Lane P's pca32 fold model (identical function, corpus and seed, nothing removed) if it
    exists, else this lane's own reference trained by the same function."""
    p = os.path.join(P_MODELS, "pca32_fold%d_s%d.pt" % (fold, seed))
    if os.path.exists(p):
        return p
    return model_path(fold, "ref", seed)


def train_without(fold, exclude_pdbs, tag, epochs=40, seed=0, verbose=True):
    """Retrain the pca32 fold model through s26/p_ladder.py's exact path with the named chains
    removed from the training corpus (peptides matched by pdb id AND sequence; fragments are
    never removed).  One .pt is the checkpoint.  Reads native distances as labels only, as the
    pinned models did (L16b)."""
    import p_ladder as PL
    path = model_path(fold, tag, seed)
    if os.path.exists(path):
        if verbose:
            print("  exists:", path)
        return path
    PL.install_guard()
    excl = {}
    for pdb in exclude_pdbs:
        rec = D.by_pdb(pdb)
        assert rec is not None, "%s is not a database peptide" % pdb
        excl[pdb] = rec.seq
    orig = PL.train_entries

    def filtered(f, n_folds=5):
        ents = orig(f, n_folds)
        kept = [e for e in ents if not (getattr(e, "pdb", None) in excl and e.seq == excl[e.pdb])]
        removed = len(ents) - len(kept)
        assert removed == len(excl), ("expected to remove %d chain(s) from fold %d's corpus, removed %d"
                                      % (len(excl), f, removed))
        return kept

    PL.train_entries = filtered
    try:
        t0 = time.time()
        data = PL.build_training("pca32", fold, verbose=verbose)
    finally:
        PL.train_entries = orig
    d_in = PL.d_in_for("pca32")
    assert data["X"].shape[1] == d_in, (data["X"].shape, d_in)
    m = PR.MLP(d_in, seed=seed, dropout=0.0, width=384, depth=3).fit(data["X"], data["Y"], epochs=epochs, verbose=verbose)
    secs = time.time() - t0
    PL.save_model(m, path, {"rung": "pca32", "fold": fold, "seed": seed, "epochs": epochs, "d_in": d_in,
                            "cfg": dict(width=384, depth=3), "n_pairs": int(len(data["Y"])),
                            "excluded": list(exclude_pdbs), "tag": tag, "train_secs": round(secs, 1),
                            "provenance": ST.provenance(__file__)})
    if verbose:
        print("  fold %d without %s: %d pairs, trained in %.0fs -> %s" % (fold, list(exclude_pdbs), len(data["Y"]), secs, path))
    del data
    return path


def train(which, epochs=40, seed=0):
    """`which`: 'ref:<fold>', 'out:<pdb>:<fold>', or 'chain' (the fixed PREREG order)."""
    if which == "chain":
        for pdb, fold in TRAIN_ORDER:
            train_without(fold, [pdb], "out_%s" % pdb, epochs=epochs, seed=seed)
        return
    kind, rest = which.split(":", 1)
    if kind == "ref":
        return train_without(int(rest), [], "ref", epochs=epochs, seed=seed)
    if kind == "out":
        pdb, fold = rest.split(":")
        return train_without(int(fold), [pdb], "out_%s" % pdb, epochs=epochs, seed=seed)
    raise ValueError(which)


def retrained_posterior(path, seq, fold):
    import p_ladder as PL
    m, meta = PL.load_model(path)
    X, i, j = PL.features_for("pca32", seq, fold)
    return m.predict_proba(X).astype(np.float64), np.asarray(i), np.asarray(j), meta


def posterior(verbose=True):
    """Native-free Part B: for each self-copy target, the reference / carrier-out / control-out
    posteriors, their differences, and the emissions with the triangle bounds."""
    tg = targets_by_pdb()
    rows = []
    for pdb, carrier, fold, _cf in SELF:
        t = tg[pdb]; seq, n = t["seq"], t["n"]
        u = load_blind(pdb)
        W = u["W"]; order = np.asarray(u["order"], int); pool = order[:K]; Wp = W[pool]
        dg = I.distogram(pdb); P0 = np.asarray(dg["prob"], float)
        e_pinned = emit(Wp, dg, seq, fold)
        arms = {"pinned": {"path": None, "prob": P0, "emit": e_pinned}}
        cands = [("reference", reference_model_path(fold)), ("carrier_out", model_path(fold, "out_%s" % carrier))]
        cands += [("control_out_%s" % c, model_path(fold, "out_%s" % c)) for c in CONTROL_CHAINS[fold]]
        for name, path in cands:
            if not os.path.exists(path):
                if verbose:
                    print("  %s: %s missing (%s); skipped" % (pdb, name, path))
                continue
            prob, i, j, meta = retrained_posterior(path, seq, fold)
            assert np.array_equal(i, np.asarray(dg["i"])) and np.array_equal(j, np.asarray(dg["j"]))
            rt = risk_from_posterior(seq, prob, i, j)
            arms[name] = {"path": os.path.relpath(path, ROOT), "prob": prob, "emit": emit(Wp, rt, seq, fold),
                          "excluded": meta.get("excluded"), "n_pairs": meta.get("n_pairs")}
        ref = arms.get("reference", arms["pinned"])
        row = {"pdb": pdb, "carrier": carrier, "fold": fold, "n": n, "reference_is": "reference" if "reference" in arms else "pinned",
               "gate": production_gate(pdb, e_pinned), "arms": {}}
        for name, a in arms.items():
            e = a["emit"]
            same_r = bool(set(e["argmin_set"].tolist()) == set(ref["emit"]["argmin_set"].tolist()))
            row["arms"][name] = {"path": a["path"], "excluded": a.get("excluded"), "n_pairs": a.get("n_pairs"),
                                 **posterior_stats(a["prob"], ref["prob"]),
                                 "vs_pinned": posterior_stats(a["prob"], P0),
                                 "top75_overlap_vs_ref": float(len(set(e["top"].tolist()) & set(ref["emit"]["top"].tolist())) / TOPM),
                                 "argmin_same_vs_ref": same_r,
                                 **triangle(ref["emit"], e, Wp, Wp, same_r),
                                 "emissions": {k: e[k] for k in ("cloud", "chain", "fit")},
                                 "argmin_windows": Wp[e["argmin_set"]]}
        assert finite(row), pdb
        rows.append(row)
        if verbose:
            print("  %s (carrier %s, fold %d) ref=%s: %s" % (pdb, carrier, fold, row["reference_is"],
                  " | ".join("%s tri_arm %.3f dE %.3f ovl %.2f" % (k, v["tri_arm"], v["expected_meanabs"], v["top75_overlap_vs_ref"])
                             for k, v in row["arms"].items())), flush=True)
        del u
    p = save("posterior", {"label": "Part B posteriors, native-free", "rows": rows}, rows=rows, n_expected=len(SELF),
             complete_keys=("pdb", "arms"))
    print("  ->", p)
    return rows


# ============================================================ GATED: endpoints
def rmsd_to_native(arrs, nat):
    """RMSD-to-native of stored emissions; the stored key `chain` is the `arm` basis."""
    return {("arm" if k == "chain" else k): float(I.ca_rmsd(np.asarray(v, float), nat)) for k, v in arrs.items()}


def sel_rmsd(windows, nat):
    """Mean RMSD-to-native over the argmin tie set (ST.argmin_tied's rule)."""
    w = np.asarray(windows, float)
    return float(np.mean([I.ca_rmsd(x, nat) for x in w]))


def endpoint(verbose=True):
    require_signoff("endpoint")
    tg = targets_by_pdb()
    out = {"label": "GATED endpoints from stored emissions", "basis": "sel 3.4540 / cloud 3.0483 / arm 3.2148 / fit 3.2041"}
    nat_of = {}

    def nat(pdb):
        if pdb not in nat_of:
            nat_of[pdb] = np.asarray(I.load_univ(pdb)["nat_ca"], float)
        return nat_of[pdb]

    # ---- Part A
    A = load_result("retrieval")
    if A and A.get("complete"):
        rowsA = []
        for r in A["rows"]:
            N = nat(r["pdb"])
            w = rmsd_to_native(r["emissions"]["with"], N); wo = rmsd_to_native(r["emissions"]["without"], N)
            w["sel"] = sel_rmsd(r["argmin_windows"]["with"], N); wo["sel"] = sel_rmsd(r["argmin_windows"]["without"], N)
            row = {"pdb": r["pdb"], "fold": r["fold"], "self_copy": r["self_copy"], "with": w, "without": wo,
                   "delta": {k: w[k] - wo[k] for k in w}, "tri": {k: r["tri_" + k] for k in ("cloud", "arm", "fit", "sel")},
                   "dropped_in_top75": r["dropped_in_top75"]}
            row["delta"]["gain"] = row["delta"]["arm"] - row["delta"]["sel"]
            if "ge06" in r:
                g = rmsd_to_native(r["ge06"]["emissions"], N); g["sel"] = sel_rmsd(r["ge06"]["argmin_windows"], N)
                row["ge06"] = {"without": g, "delta": {k: w[k] - g[k] for k in g}, "n_dropped": r["ge06"]["n_dropped"]}
            ins = []
            for x in r["oracle_insertion"]:
                e = rmsd_to_native(x["emissions"], N); e["sel"] = sel_rmsd(x["argmin_windows"], N)
                ins.append({"carrier_len": x["carrier_len"], "inserted_in_top75": x["inserted_in_top75"],
                            "inserted_window_rr": float(I.ca_rmsd(np.asarray(x["inserted_window"], float), N)),
                            "with_insertion": e, "delta_insert_minus_prod": {k: e[k] - w[k] for k in e}})
            row["oracle_insertion"] = ins
            rowsA.append(row)
        selfA = [r for r in rowsA if r["self_copy"]]; ctrlA = [r for r in rowsA if not r["self_copy"]]
        for b in ("sel", "cloud", "arm", "fit", "gain"):
            for r in selfA:
                cd = np.abs([c["delta"][b] for c in ctrlA]); r.setdefault("ctrl_pctile", {})[b] = float((cd <= abs(r["delta"][b])).mean())
        out["A"] = {"rows": rowsA,
                    "self": {r["pdb"]: {"delta": r["delta"], "tri": r["tri"], "ctrl_pctile": r["ctrl_pctile"]} for r in selfA},
                    "control_abs_delta_pctiles": {b: {q: float(np.percentile(np.abs([c["delta"][b] for c in ctrlA]), q)) for q in (50, 90, 100)}
                                                  for b in ("sel", "cloud", "arm", "fit", "gain")},
                    "control_mean_delta": {b: float(np.mean([c["delta"][b] for c in ctrlA])) for b in ("sel", "cloud", "arm", "fit")},
                    "ge06_mean_delta_all126": {b: float(np.mean([r["ge06"]["delta"][b] if "ge06" in r else 0.0 for r in rowsA]))
                                               for b in ("sel", "cloud", "arm", "fit")},
                    "ge06_n_targets": int(sum(1 for r in rowsA if "ge06" in r)),
                    "insertion": [{"pdb": r["pdb"], **x} for r in rowsA for x in r["oracle_insertion"]]}
        if verbose:
            print("PART A (channel A), signed delta = with - without (positive: the leaked emission is WORSE than the un-leaked one)")
            for r in selfA:
                print("  %s  delta arm %+.4f cloud %+.4f sel %+.4f fit %+.4f gain %+.4f | tri arm %.3f | ctrl pctile arm %.2f"
                      % (r["pdb"], r["delta"]["arm"], r["delta"]["cloud"], r["delta"]["sel"], r["delta"]["fit"], r["delta"]["gain"],
                         r["tri"]["arm"], r["ctrl_pctile"]["arm"]))
            print("  ge06 (S10-4 operator) mean delta over 126 on arm %+.5f, fit %+.5f, sel %+.5f (n targets %d)"
                  % (out["A"]["ge06_mean_delta_all126"]["arm"], out["A"]["ge06_mean_delta_all126"]["fit"],
                     out["A"]["ge06_mean_delta_all126"]["sel"], out["A"]["ge06_n_targets"]))
    # ---- Part C
    C = load_result("envelope")
    if C and C.get("complete"):
        rowsC = []
        for r in C["rows"]:
            N = nat(r["pdb"])
            clean = rmsd_to_native(r["clean"], N); clean["sel"] = sel_rmsd(r["clean_argmin_windows"], N)
            leaked = {}
            for j, v in r["leaked"].items():
                e = rmsd_to_native(v["emissions"], N); e["sel"] = sel_rmsd(v["argmin_windows"], N)
                leaked[j] = e
            rowsC.append({"pdb": r["pdb"], "fold": r["fold"], "clean": clean, "leaked": leaked,
                          "leaked_mean": {b: float(np.mean([leaked[j][b] for j in leaked])) for b in clean},
                          "leaked_spread": {b: float(np.std([leaked[j][b] for j in leaked])) for b in clean}})
        pdbs = [r["pdb"] for r in rowsC]; folds = ST.pinned_folds(pdbs)
        stats = {}
        for b in ("arm", "cloud", "sel", "fit"):
            a = np.array([r["leaked_mean"][b] for r in rowsC]); c = np.array([r["clean"][b] for r in rowsC])
            res = ST.compare(a, c, folds, names=pdbs, label="Part C envelope: leaked-model mean minus clean [%s] (ORACLE)" % b)
            stats[b] = res
            if verbose:
                print(ST.fmt(res))
            per = {}
            for j in range(5):
                idx = [k for k, r in enumerate(rowsC) if str(j) in r["leaked"]]
                aj = np.array([rowsC[k]["leaked"][str(j)][b] for k in idx]); cj = c[idx]
                per[j] = ST.compare(aj, cj, folds[idx], names=[pdbs[k] for k in idx], label="  model %d minus clean [%s]" % (j, b))
            stats[b + "_per_model"] = per
        out["C"] = {"rows": rowsC, "stats": stats,
                    "leaked_spread_mean": {b: float(np.mean([r["leaked_spread"][b] for r in rowsC])) for b in ("arm", "cloud", "sel", "fit")}}
    # ---- Part B
    B = load_result("posterior")
    if B and B.get("complete"):
        rowsB = []
        for r in B["rows"]:
            N = nat(r["pdb"])
            arms = {}
            for name, a in r["arms"].items():
                e = rmsd_to_native(a["emissions"], N); e["sel"] = sel_rmsd(a["argmin_windows"], N); e["gain"] = e["arm"] - e["sel"]
                arms[name] = e
            ref = arms.get("reference", arms["pinned"])
            rowsB.append({"pdb": r["pdb"], "carrier": r["carrier"], "fold": r["fold"], "reference_is": r["reference_is"], "arms": arms,
                          "delta_vs_ref": {name: {b: ref[b] - arms[name][b] for b in ref} for name in arms}})
        out["B"] = {"rows": rowsB}
        if verbose:
            print("PART B (channel B), delta = reference - arm (positive: removing that chain made the answer WORSE, i.e. the chain helped)")
            for r in rowsB:
                for name, d in r["delta_vs_ref"].items():
                    print("  %s %-22s arm %+.4f cloud %+.4f sel %+.4f gain %+.4f" % (r["pdb"], name, d["arm"], d["cloud"], d["sel"], d["gain"]))
    # ---- both channels removed on the 4 (the clean counterfactual): carrier-out model x self-window dropped
    if A and A.get("complete") and B and B.get("complete"):
        both = []
        for rB in B["rows"]:
            pdb = rB["pdb"]; t = tg[pdb]; N = nat(pdb)
            if "carrier_out" not in rB["arms"]:
                continue
            u = I.load_univ(pdb); W = u["W"]; order = np.asarray(u["order"], int)
            exact = self_windows(u, t["seq"]); pool1 = drop_and_refill(order, exact)
            path = model_path(t["fold"], "out_%s" % rB["carrier"])
            prob, i, j, _ = retrained_posterior(path, t["seq"], t["fold"])
            rt = risk_from_posterior(t["seq"], prob, i, j)
            e = emit(W[pool1], rt, t["seq"], t["fold"])
            clean = rmsd_to_native({k: e[k] for k in ("cloud", "chain", "fit")}, N); clean["sel"] = sel_rmsd(W[pool1][e["argmin_set"]], N)
            clean["gain"] = clean["arm"] - clean["sel"]
            prodA = dict([x for x in out["A"]["rows"] if x["pdb"] == pdb][0]["with"]); prodA["gain"] = prodA["arm"] - prodA["sel"]
            both.append({"pdb": pdb, "production": prodA, "both_removed": clean,
                         "delta": {b: prodA[b] - clean[b] for b in ("sel", "cloud", "arm", "fit", "gain")}})
        out["both_removed"] = both
        if verbose:
            for r in both:
                print("  BOTH REMOVED %s: delta arm %+.4f cloud %+.4f sel %+.4f gain %+.4f" % (r["pdb"], r["delta"]["arm"], r["delta"]["cloud"], r["delta"]["sel"], r["delta"]["gain"]))
    p = save("endpoint", out)
    print("  ->", p)
    return out


def floor(verbose=True):
    """GATED Part E: the same sequence in a different deposit, the 18 verbatim dev cases."""
    require_signoff("floor")
    tg = I.targets(); folds = D.folds(5); dev = {t["pdb"] for t in tg}
    rows = []
    for t in tg:
        seq, n = t["seq"], t["n"]
        rec = D.by_pdb(t["pdb"]); assert rec is not None and rec.seq == seq
        parts = [p for p in D.load() if p.seq != seq and (seq in p.seq or p.seq in seq)]
        if not parts:
            continue
        u = I.load_univ(t["pdb"]); N = np.asarray(u["nat_ca"], float); rr = np.asarray(u["rr"], float)[np.asarray(u["order"], int)[:K]]
        db_matches_instrument = bool(np.allclose(np.asarray(rec.ca, float), N))   # recorded, not asserted
        for p in parts:
            if seq in p.seq:
                s = p.seq.find(seq); a, b = N, np.asarray(p.ca[s:s + n], float); role = "carrier"
            else:
                s = seq.find(p.seq); a, b = N[s:s + len(p.seq)], np.asarray(p.ca, float); role = "carried"
            rows.append({"pdb": t["pdb"], "n": n, "fold": t["fold"], "partner": p.pdb if p.pdb in dev else "(non-dev n=%d)" % len(p.seq),
                         "partner_len": len(p.seq), "role": role, "start": int(s), "same_fold": folds[p.seq] == t["fold"],
                         "longer_identity": float(D.identity(seq, p.seq)), "overlap_len": int(min(n, len(p.seq))),
                         "db_ca_matches_instrument_native": db_matches_instrument,
                         "cross_deposit_rmsd": float(I.ca_rmsd(a, b)),
                         "pool_mean_rr": float(rr.mean()), "pool_best_rr": float(rr.min())})
    vals = np.array([r["cross_deposit_rmsd"] for r in rows])
    summ = {"n_pairs": len(rows), "n_targets": len({r["pdb"] for r in rows}), "median": float(np.median(vals)), "mean": float(vals.mean()),
            "min": float(vals.min()), "max": float(vals.max()), "frac_below_1.0": float((vals < 1.0).mean()), "frac_below_1.5": float((vals < 1.5).mean()),
            "cross_fold": {r["pdb"]: r["cross_deposit_rmsd"] for r in rows if not r["same_fold"]},
            "pool_mean_rr_mean": float(np.mean([r["pool_mean_rr"] for r in rows])), "ensemble_spread_record": 1.044}
    if verbose:
        for r in rows:
            print("  %s n=%2d %s %-16s len %2d %s id %.3f  cross-deposit %.3f  (pool mean %.3f best %.3f)"
                  % (r["pdb"], r["n"], r["role"], r["partner"], r["partner_len"], "SAME " if r["same_fold"] else "OTHER", r["longer_identity"],
                     r["cross_deposit_rmsd"], r["pool_mean_rr"], r["pool_best_rr"]))
        print(json.dumps(summ, indent=1))
    p = save("floor", {"label": "Part E: same sequence, different deposit (ORACLE DIAGNOSTIC)", "summary": summ, "rows": rows},
             rows=rows, n_expected=len(rows), complete_keys=("pdb", "cross_deposit_rmsd"))
    print("  ->", p)
    return rows


# ============================================================ Part D: the bound
def bound_from_deltas(deltas, n_leaked=N_LEAKED_BENCH, n_bench=N_BENCH):
    d = np.abs(np.asarray(deltas, float))
    return float(n_leaked / n_bench * d.max()) if len(d) else float("nan")


def materiality(x):
    if not np.isfinite(x):
        return "not measured"
    if x < IMMATERIAL:
        return "IMMATERIAL (< %.3f A, one tenth of the benchmark CI half-width %.3f)" % (IMMATERIAL, BENCH_CI_HALF)
    if x < BENCH_CI_HALF:
        return "MINOR (%.3f to %.3f A)" % (IMMATERIAL, BENCH_CI_HALF)
    return "MATERIAL (>= %.3f A)" % BENCH_CI_HALF


def report(verbose=True):
    """Part D.  Native-free bounds (triangle) are always available once `retrieval` /
    `envelope` / `posterior` exist; the signed bounds need `endpoint` (gated)."""
    out = {"label": "the 2/60 bound", "n_bench": N_BENCH, "n_leaked_bench": N_LEAKED_BENCH, "bench_ci_half": BENCH_CI_HALF,
           "immaterial_threshold": IMMATERIAL, "assumptions": "s26/PREREG_selfcopy_bound.md section 2 (A1-A5)"}
    A = load_result("retrieval"); C = load_result("envelope"); B = load_result("posterior"); E = load_result("endpoint")
    tri = {}
    if A and A.get("complete"):
        s = [r for r in A["rows"] if r["self_copy"]]
        tri["A_real4"] = {b: bound_from_deltas([r["tri_" + b] for r in s]) for b in ("cloud", "arm", "fit", "sel")}
        tri["A_real4"]["gain"] = bound_from_deltas([r["tri_arm"] + r["tri_sel"] for r in s])
        tri["A_real4_per_target"] = {r["pdb"]: {b: r["tri_" + b] for b in ("cloud", "arm", "fit", "sel")} for r in s}
    if B and B.get("complete"):
        tri["B_carrier_out"] = {b: bound_from_deltas([r["arms"]["carrier_out"]["tri_" + b] for r in B["rows"] if "carrier_out" in r["arms"]])
                                for b in ("cloud", "arm", "fit", "sel")}
    if C and C.get("complete"):
        tri["C_envelope_max_over_models"] = {b: bound_from_deltas([max(v["tri_" + b] for v in r["leaked"].values()) for r in C["rows"]])
                                             for b in ("cloud", "arm", "fit", "sel")}
        tri["C_envelope_p95_over_targets"] = {b: float(N_LEAKED_BENCH / N_BENCH * np.percentile(
            [np.mean([v["tri_" + b] for v in r["leaked"].values()]) for r in C["rows"]], 95)) for b in ("cloud", "arm", "fit", "sel")}
    out["triangle_bounds_native_free"] = tri
    signed = {}
    if E:
        if "A" in E:
            s = E["A"]["self"]
            signed["A_real4"] = {b: bound_from_deltas([s[p]["delta"][b] for p in s]) for b in ("sel", "cloud", "arm", "fit", "gain")}
        if "both_removed" in E and E["both_removed"]:
            signed["both_removed4"] = {b: bound_from_deltas([r["delta"][b] for r in E["both_removed"]]) for b in ("sel", "cloud", "arm", "fit", "gain")}
        if "C" in E:
            st = E["C"]["stats"]
            signed["C_envelope_fold_ci"] = {b: {"effect": st[b]["effect"], "mde": st[b]["mde"], "ci95_fold": st[b]["ci95_fold"],
                                                "verdict": st[b]["verdict"],
                                                "bound": float(N_LEAKED_BENCH / N_BENCH * max(abs(x) for x in st[b]["ci95_fold"]))}
                                            for b in ("arm", "cloud", "sel", "fit")}
    out["signed_bounds_gated"] = signed
    verdict = {}
    for b in ("arm", "gain", "sel", "cloud"):
        cands = []
        if "both_removed4" in signed and b in signed["both_removed4"]:
            cands.append(signed["both_removed4"][b])
        elif "A_real4" in signed and b in signed["A_real4"]:
            cands.append(signed["A_real4"][b])
        if "C_envelope_fold_ci" in signed and b in signed["C_envelope_fold_ci"]:
            cands.append(signed["C_envelope_fold_ci"][b]["bound"])
        if not cands and "A_real4" in tri and b in tri["A_real4"]:
            cands.append(tri["A_real4"][b])
        x = max(cands) if cands else float("nan")
        verdict[b] = {"bound_A": x, "class": materiality(x), "source": "signed" if signed else "triangle (native-free)"}
    out["verdict"] = verdict
    if verbose:
        print(json.dumps({"triangle": tri, "signed": signed, "verdict": verdict}, indent=1))
    p = save("bound", out)
    print("  ->", p)
    return out


# ============================================================ CLI
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("census", "retrieval", "envelope", "train", "posterior", "endpoint", "floor", "report"))
    ap.add_argument("--probe", default=None, help="one target (pdb) for the memory probe")
    ap.add_argument("--which", default="chain", help="train: 'chain' | 'ref:<fold>' | 'out:<pdb>:<fold>'")
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    if a.cmd == "census":
        census()
    elif a.cmd == "retrieval":
        retrieval(probe=a.probe)
    elif a.cmd == "envelope":
        envelope(probe=a.probe)
    elif a.cmd == "train":
        train(a.which, epochs=a.epochs, seed=a.seed)
    elif a.cmd == "posterior":
        posterior()
    elif a.cmd == "endpoint":
        endpoint()
    elif a.cmd == "floor":
        floor()
    elif a.cmd == "report":
        report()
    return 0


if __name__ == "__main__":
    sys.exit(main())
