"""Run a set of retrieval-key ARMS end to end through the full production chain.

Each arm supplies a per-target score over the whole universe; the top-K=500 by that score
is the pool; the pool goes through shipped-score top-75 -> coordinate average -> project
(lam=0 'fit' and lam=0.3 'proj').  One pass over the 126 targets, all arms per target, so
each universe is loaded once.

Usage:  python -m s12.key_emit --arms blosum,oracle_native --out key_emit_oracle
"""
from __future__ import annotations
import os, sys, json, time, argparse
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import key_lib as KL


# ------------------------------------------------------------------ arm score functions
def make_arms(names, alph="abego4", pred=None, wmix=None):
    """Return {name: fn(u) -> score over universe}.  `pred` is {pdb: (n,nb) posterior}."""
    nat = KL.native_torsions()
    arms = {}

    def blosum(u):
        return u["sim"].astype(float)

    def _wb(u):
        return KL.bins_of(alph, u["PHI"], u["PSI"])

    def oracle_native(u):
        nb = KL.bins_of(alph, *nat[u["pdb"]])
        return KL.hard_key_score(_wb(u), nb)

    def oracle_bestwin(u):
        wb = _wb(u)
        return KL.hard_key_score(wb, wb[int(np.argmin(u["rr"]))])

    def _pred_score(u, tag):
        P = pred[tag][u["pdb"]]
        return KL.soft_key_score(_wb(u), P)

    def rand_marginal(u):
        """CAPACITY NULL: key drawn i.i.d. from the library's pooled bin marginal."""
        from s12.key_noise import MARGINAL
        nb = KL.nbins(alph)
        rng = np.random.default_rng(abs(hash(u["pdb"])) % (2 ** 31))
        p = MARGINAL[:nb] / MARGINAL[:nb].sum()
        return KL.hard_key_score(_wb(u), rng.choice(nb, u["n"], p=p))

    def oracle_ss3(u):
        nb = KL.bins_of("ss3", *nat[u["pdb"]])
        return KL.hard_key_score(KL.bins_of("ss3", u["PHI"], u["PSI"]), nb)

    reg = {"blosum": blosum, "oracle_native": oracle_native, "oracle_bestwin": oracle_bestwin,
           "rand_marginal": rand_marginal, "oracle_ss3": oracle_ss3}
    for nm in names:
        if nm in reg:
            arms[nm] = reg[nm]
        elif nm.startswith("pred:"):
            tag = nm.split(":", 1)[1]
            arms[nm] = (lambda t: (lambda u: _pred_score(u, t)))(tag)
        elif nm.startswith("mix:"):
            _, tag, w = nm.split(":")
            w = float(w)
            arms[nm] = (lambda t, w: (lambda u: (1 - w) * KL.zscore(u["sim"].astype(float))
                                      + w * KL.zscore(_pred_score(u, t))))(tag, w)
        elif nm.startswith("mixoracle:"):
            w = float(nm.split(":")[1])
            arms[nm] = (lambda w: (lambda u: (1 - w) * KL.zscore(u["sim"].astype(float))
                                  + w * KL.zscore(oracle_native(u))))(w)
        else:
            raise ValueError(nm)
    return arms


def run(arm_names, out, alph="abego4", pred_path=None, pdbs=None, lam=0.3, m=I.M):
    pred = {}
    if pred_path:
        with open(pred_path) as fh:
            raw = json.load(fh)
        for tag, d in raw.items():
            pred[tag] = {k: np.asarray(v, float) for k, v in d.items()}
    arms = make_arms(arm_names, alph=alph, pred=pred)
    tg = {t["pdb"]: t for t in I.targets()}
    pdbs = pdbs or list(tg)
    rows = []
    outp = os.path.join(I.RESULTS, out + ".json")
    t0 = time.time()
    for c, pdb in enumerate(pdbs):
        while I.free_gb() < 1.2:
            time.sleep(20)
        u = I.load_univ(pdb); dg = I.distogram(pdb)
        rr = u["rr"]
        row = {"pdb": pdb, "n": u["n"], "fold": u["fold"], "f18": pdb in I.FAIL18,
               "univ_best": float(rr.min()), "arms": {}}
        for nm, fn in arms.items():
            idx = KL.retrieve(fn(u))
            st = KL.pool_stats(rr, idx)
            e = KL.emit(u, idx, dg, m=m, lam=lam)
            row["arms"][nm] = {**st, **e}
        rows.append(row)
        del u, dg
        with open(outp, "w") as fh:
            json.dump(rows, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
        print(f"[{c+1}/{len(pdbs)}] {pdb} " +
              " ".join(f"{k}:fit={v['fit']:.3f}/best={v['best']:.3f}" for k, v in row["arms"].items()) +
              f"  t={time.time()-t0:.0f}s free={I.free_gb():.2f}", flush=True)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--alph", default="abego4")
    ap.add_argument("--pred", default=None)
    ap.add_argument("--n", type=int, default=0)
    a = ap.parse_args()
    tg = [t["pdb"] for t in I.targets()]
    run(a.arms.split(","), a.out, alph=a.alph, pred_path=a.pred, pdbs=tg[: a.n] if a.n else tg)
