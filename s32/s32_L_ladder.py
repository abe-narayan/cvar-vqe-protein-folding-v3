"""LANE L / L2 -- the loss ladder at 40-60 residues, and the SAME ladder at 9-16.

PREREG `s32/PREREG_S32_L.md` at commit 88f2da39, hypothesis L-H2, predictions P1-P3.

THE LIKE-FOR-LIKE PROBLEM, and how it is solved
The canonical ladder's middle rungs are defined by the DISTOGRAM score's top-128 and
top-75.  That score cannot be evaluated at 40-60 residues, and the reason is structural
rather than a matter of accuracy: `core/predict.py` hard-codes `MAXLEN = 26`, its
separation one-hot `SEP_BINS` tops out at 24 so every pair beyond that collapses into a
single bin, and both `sep / 26.0` and `n / 26.0` are normalisations by 26.  At n = 55 the
raw `n` feature and 27% of all pairs sit outside the range the model was ever fitted on.

So this module does NOT compare a new ladder against the old one.  It builds **one ladder
that is length-portable at every rung** -- BLOSUM62 retrieval, ORACLE best member, ORACLE
sparse convex combination, uniform averaging -- and runs that SAME ladder, through the SAME
code, at BOTH lengths.  The comparison is then length against length, with the operator
held fixed, which is the only comparison that can answer whether the ladder's shape is a
peptide artefact (contract rules 4 and 6).

RUNGS, all measured on the BUILT CHAIN (the projector's output).  Basis is stated on every
number: **built chain, pre-AMBER**.  Stage 4 is a validity stage priced by S8-12 at +0.011
to +0.026 A on the canonical instrument; it is omitted at BOTH lengths, so every rung
DIFFERENCE reported here is unaffected by its absence.

    A  pool_best      ORACLE / NOT DEPLOYABLE   best single member of the K=500 pool
    B  sparse_s10     ORACLE / NOT DEPLOYABLE   best sparse convex combination, s = 10
    C  top75_best     ORACLE / NOT DEPLOYABLE   best single member of the BLOSUM top-75
    D  avg75          DEPLOYABLE-SHAPED         uniform average of the BLOSUM top-75
    E  avg75_random   ZERO-INFORMATION CONTROL  uniform average of 75 members drawn at
                                                random from the K=500 pool

E is the control for D and it is matched to D's own space (contract rule 6): the same
operator, the same member count, the same pool, differing only in whether the 75 were
chosen by the retrieval key.  It is drawn `N_DRAW` times per target and the DRAW MEAN and
draw-to-draw sd are both reported, never the best draw (contract rule 10).

    python -m s32.s32_L_ladder long  --workers 2
    python -m s32.s32_L_ladder short --workers 2
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "s32", "results")
os.makedirs(RESULTS, exist_ok=True)

K = 500
M = 75
SPARSE_S = 10
N_DRAW = 3
RNG_SEED = 0


# ----------------------------------------------------------------------------- geometry
def _ins():
    sys.path.insert(0, BASE)
    from s12 import instrument as I
    return I


def _project(C, seq, fold):
    """STAGE 3b exactly as production: lambda path to 0.3, multi-start, exact gradient."""
    from core import project as pj
    pen = pj.make_penalty("ramah", seq, int(fold))
    path = pj.lam_path(np.asarray(C, float), pen, (0.0, 0.3), maxiter=300,
                       multi=True, grad="exact")
    return np.asarray(path[0.3][0], float)


def _avg(W, idx):
    I = _ins()
    Ws = np.asarray(W, float)[np.asarray(idx, int)]
    P = I.pairwise_rmsd(Ws)
    C, _ = I.coordinate_average(Ws, P)
    return C


# ------------------------------------------------------------------------ the candidate
def _bank_long(t, corpus=None, exclude_kmer=9):
    """Top-K BLOSUM62 length-n windows of the cached corpus, leakage-filtered.

    The corpus is loaded from the module-level cache INSIDE the worker.  Passing it as a
    task argument would pickle ~100 MB per target through the pool and is how a run of this
    shape exhausts a box that has 2.7 GB free.

    LEAKAGE.  The first implementation ran a Needleman-Wunsch identity of the target against
    every one of 13,751 chains.  That is a 55 x 2000 dynamic program per chain in Python and
    it cost 600 s per target, which is 7.6 hours for the instrument.  It is also the WRONG
    statistic: `s32_L_leaknull` measured that alignment identity against short bank members
    is at the null here, while the verbatim-substring statistic separates real from shuffled
    perfectly (real 0.330, shuffled 0.000).  So the filter is the one with the measured null:
    **a chain is dropped if it shares any verbatim 9-mer with the target.**  Homologs share
    many; unrelated chains share none.  It is both correct and ~300x cheaper.

    SPEED.  Windows are cut with `sliding_window_view` over pre-encoded chains and scored by
    one gather-and-sum per chain, instead of one `encode` call per window.  Only window
    REFERENCES (chain, start, score) are accumulated; coordinates are materialised for the
    top K alone.
    """
    import s7.audit as audit
    if corpus is None:
        from s32 import s32_L_corpus as cp
        corpus = cp.load()
    n = int(t["n"]); tseq = t["seq"]; code = audit.encode(tseq)
    tk = {tseq[i:i + exclude_kmer] for i in range(len(tseq) - exclude_kmer + 1)}
    own = t["pdb"].upper()
    swv = np.lib.stride_tricks.sliding_window_view

    ci, st, sims = [], [], []
    for a in range(len(corpus["pdb"])):
        if corpus["pdb"][a] == own:
            continue
        seq = corpus["seq"][a]
        m = len(seq)
        if m <= n + 1:
            continue
        if any(k in seq for k in tk):            # homolog (or the target itself) -- drop
            continue
        ok = corpus["ok"][a]                     # CA step in (3.5, 4.1) for bond s -> s+1
        fin = corpus["fin"][a]                   # finite phi/psi and a real residue code
        cok = np.concatenate([[0.0], np.cumsum(~ok)])
        cfin = np.concatenate([[0.0], np.cumsum(~fin)])
        s0 = np.arange(1, m - n)
        # a window [s, s+n) is usable iff no broken bond inside it and no bad residue in it
        good = ((cok[s0 + n - 1] - cok[s0]) == 0) & ((cfin[s0 + n] - cfin[s0]) == 0)
        if not good.any():
            continue
        C = corpus["code"][a]
        Wc = swv(C, n)[s0[good]]                 # (nw, n) int codes, a VIEW
        sm = audit.B62[Wc, code[None, :]].sum(1)
        ci.append(np.full(len(sm), a, np.int32))
        st.append(s0[good].astype(np.int32))
        sims.append(sm.astype(np.float32))
    ci = np.concatenate(ci); st = np.concatenate(st); sims = np.concatenate(sims)
    order = np.argsort(-sims, kind="stable")[:K]
    W = np.stack([np.asarray(corpus["ca"][ci[i]], float)[st[i]:st[i] + n] for i in order])
    src = [f"{corpus['pdb'][ci[i]]}_{int(st[i])}" for i in order]
    return W, sims[order].astype(float), src


# ----------------------------------------------------------------------------- one target
def run_target(kind, t, corpus=None):
    I = _ins()
    t0 = time.time()
    if kind == "long":
        W, sim, src = _bank_long(t, corpus)
        from core import geometry as geo
        import glob
        p = os.path.join(BASE, "prots", t["pdb"] + ".pdb")
        if not os.path.exists(p):
            p = glob.glob(os.path.join(BASE, "prots", t["pdb"].lower() + ".pdb"))[0]
        seq, coords, _, _ = geo.native_coords_from_pdb(p)
        nat = np.asarray(coords["CA"], float)
        fold = int(t["fold"])
    else:
        u = I.load_univ(t["pdb"])
        pool = I.pool_idx(u, K)
        W = np.asarray(u["W"], float)[pool]
        sim = np.asarray(u["sim"], float)[pool]
        nat = np.asarray(u["nat_ca"], float)
        seq = t["seq"]; fold = int(t["fold"])
    n = len(nat)
    rr = I.kabsch_rmsd_batch(W, nat)                                        # ORACLE labels
    order = np.arange(len(W))                       # W is already in BLOSUM (stable) order

    out = {"pdb": t["pdb"], "n": n, "fold": fold, "kind": kind, "K": len(W),
           "cloud": {}, "chain": {}}

    # A -- ORACLE best single member of the pool
    a = int(np.argmin(rr))
    out["cloud"]["pool_best"] = float(rr[a])
    out["chain"]["pool_best"] = I.ca_rmsd(_project(W[a], seq, fold), nat)

    # C -- ORACLE best single member of the BLOSUM top-75
    c = int(np.argmin(rr[:M]))
    out["cloud"]["top75_best"] = float(rr[c])
    out["chain"]["top75_best"] = I.ca_rmsd(_project(W[c], seq, fold), nat)

    # B -- ORACLE best sparse convex combination, s = 10, over the K=500 pool
    from s29.s29_O_ladder import oracle_sparse_greedy, superpose_one
    Wp = np.stack([superpose_one(w, W[0]) for w in W]).reshape(len(W), -1)
    sp = oracle_sparse_greedy(Wp, nat, n, s_list=(SPARSE_S,))
    Xs = sp[SPARSE_S][1]
    out["cloud"]["sparse_s10"] = float(sp[SPARSE_S][0])
    out["chain"]["sparse_s10"] = I.ca_rmsd(_project(Xs, seq, fold), nat)

    # D -- the deployable-shaped readout: uniform average of the BLOSUM top-75
    Cd = _avg(W, order[:M])
    out["cloud"]["avg75"] = I.ca_rmsd(Cd, nat)
    out["chain"]["avg75"] = I.ca_rmsd(_project(Cd, seq, fold), nat)

    # E -- ZERO-INFORMATION CONTROL, matched to D's own space; the DRAW MEAN is reported
    rng = np.random.default_rng(RNG_SEED + abs(hash(t["pdb"])) % 10_000)
    dc, dh = [], []
    for _ in range(N_DRAW):
        idx = rng.choice(len(W), M, replace=False)
        Ce = _avg(W, idx)
        dc.append(I.ca_rmsd(Ce, nat))
        dh.append(I.ca_rmsd(_project(Ce, seq, fold), nat))
    out["cloud"]["avg75_random"] = float(np.mean(dc))
    out["chain"]["avg75_random"] = float(np.mean(dh))
    out["draw_sd"] = {"cloud": float(np.std(dc, ddof=1)), "chain": float(np.std(dh, ddof=1))}
    out["draws"] = {"cloud": dc, "chain": dh}
    out["secs"] = time.time() - t0
    return out


# ---------------------------------------------------------------------------- the driver
def _targets(kind, limit=None):
    if kind == "long":
        tg = json.load(open(os.path.join(RESULTS, "long40_manifest.json")))["targets"]
    else:
        tg = _ins().targets()
    return tg[:limit] if limit else tg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["long", "short"])
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(v, "1")
    tg = _targets(a.kind, a.limit or None)
    outp = os.path.join(RESULTS, f"L2_ladder_{a.kind}.jsonl")
    done = set()
    if os.path.exists(outp):
        for line in open(outp):
            try:
                done.add(json.loads(line)["pdb"])
            except Exception:
                pass
    todo = [t for t in tg if t["pdb"] not in done]
    print(f"{a.kind}: {len(todo)}/{len(tg)} targets to run, {a.workers} workers", flush=True)
    corpus = None      # loaded inside each worker; see `_bank_long`
    t0 = time.time()
    if a.workers <= 1:
        for k, t in enumerate(todo):
            r = run_target(a.kind, t, corpus)
            with open(outp, "a") as fh:
                fh.write(json.dumps(r) + "\n")
            print(f"[{k+1}/{len(todo)}] {t['pdb']} n={r['n']} "
                  f"pool_best={r['chain']['pool_best']:.3f} "
                  f"avg75={r['chain']['avg75']:.3f} ({r['secs']:.0f}s)", flush=True)
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            futs = {ex.submit(run_target, a.kind, t, corpus): t for t in todo}
            for k, f in enumerate(as_completed(futs)):
                try:
                    r = f.result()
                except Exception as e:                                   # noqa: BLE001
                    print(f"  FAIL {futs[f]['pdb']}: {e}", flush=True)
                    continue
                with open(outp, "a") as fh:
                    fh.write(json.dumps(r) + "\n")
                print(f"[{k+1}/{len(todo)}] {r['pdb']} n={r['n']} "
                      f"pool_best={r['chain']['pool_best']:.3f} "
                      f"avg75={r['chain']['avg75']:.3f} ({r['secs']:.0f}s)", flush=True)
    print(f"done in {time.time()-t0:.0f}s -> {outp}", flush=True)


if __name__ == "__main__":
    main()
