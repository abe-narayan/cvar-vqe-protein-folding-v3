"""SPRINT 15, INFO, PART B.3 -- does the in-band sign result TRANSFER to the real pools?

`s15.info_chain` found the first demonstrated native-free per-target sign signal: inside
the near-native band of the enumerated k=4 space, the distogram-predicted radius of
gyration picks the sign of the compactness ordering on 18 of 19 targets and captures 99.7%
of the ORACLE sign channel.

The project's own ledger says a decoy bank is not a pool proxy and that skill measured on a
constructed space has repeatedly failed to transfer.  So the claim is re-run, unchanged in
logic, on the actual 126-target retrieval instrument:

    space     the shipped K=500 BLOSUM window pool for each of the 126 tuning targets
    band      pool members within `pool_best + 1.5 A` (the instrument's own BAND)
    objective E_s = s * Rg(window),  s in {-1, +1}
    operator  top-m mean CA-RMSD, m = 24 and 75 (the shipped sub-selection is 75)
    rule      the sign is chosen LEAVE-ONE-TARGET-OUT from a native-free proxy; the
              threshold and the orientation are both fitted on the other 125 targets

Nulls: a permuted-proxy arm (labels shuffled across targets), the majority/fixed-sign rule,
the exact coin-flip average, and the ORACLE sign ceiling.

    python -m s15.info_signpool
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I          # noqa: E402
from s15 import info_lib as L           # noqa: E402
from s15.info_chain import emit          # noqa: E402

BAND = 1.5


def build(cache=os.path.join(L.CACHE, "signpool.json")):
    if os.path.exists(cache):
        with open(cache) as fh:
            return json.load(fh)
    with open(os.path.join(ROOT, "s14", "results", "signpred.json")) as fh:
        sp = {r["pdb"]: r for r in json.load(fh)["per_target"]}
    rows = []
    for t in I.targets():
        p, n, fold = t["pdb"], int(t["n"]), int(t["fold"])
        u = I.load_univ(p)
        pi = I.pool_idx(u)
        W = np.asarray(u["W"], float)[pi]
        rr = np.asarray(u["rr"], float)[pi]
        rg = L.rg_of(W)
        band = np.flatnonzero(rr <= rr.min() + BAND)
        r = {"pdb": p, "n": n, "fold": fold,
             "pool_mean": float(rr.mean()), "pool_best": float(rr.min()),
             "band_n": int(len(band)), "band_mean": float(rr[band].mean()),
             "band_best": float(rr[band].min()),
             "rg_pool_mean": float(rg.mean()), "rg_pool_sd": float(rg.std()),
             "rg_pool_iqr": float(np.percentile(rg, 75) - np.percentile(rg, 25)),
             "rg_nat_ORACLE": L.rg_of(np.asarray(u["nat_ca"], float)[None])[0]}
        for m in (24, 75):
            r[f"p{m}_plus"] = emit(rg, rr, m)
            r[f"p{m}_minus"] = emit(-rg, rr, m)
            mm = min(m, len(band))
            r[f"b{m}_plus"] = emit(rg[band], rr[band], mm)
            r[f"b{m}_minus"] = emit(-rg[band], rr[band], mm)
        s = sp[p]
        for k in ("rg_disto", "rg_pool", "rg_emit"):
            r[k] = float(s[k])
            r["z_" + k] = float((s[k] - rg.mean()) / rg.std())
        r["z_rg_nat_ORACLE"] = float((r["rg_nat_ORACLE"] - rg.mean()) / rg.std())
        rows.append(r)
    with open(cache, "w") as fh:
        json.dump(rows, fh)
    return rows


def loo_sign_rule(x, better_minus):
    x = np.asarray(x, float); bm = np.asarray(better_minus, bool)
    n = len(x)
    out = np.empty(n, int)
    for i in range(n):
        k = np.ones(n, bool); k[i] = False
        thr = np.median(x[k])
        hi = x[k] > thr
        agree = (hi == bm[k]).mean()
        ph = x[i] > thr
        out[i] = (-1 if ph else +1) if agree >= 0.5 else (+1 if ph else -1)
    return out


def loo_majority(better_minus):
    bm = np.asarray(better_minus, bool)
    n = len(bm)
    out = np.empty(n, int)
    for i in range(n):
        k = np.ones(n, bool); k[i] = False
        out[i] = -1 if bm[k].mean() > 0.5 else +1
    return out


def axis(rows, pk, mk, name, base_key, n_perm=2000, seed=0):
    pdbs = [r["pdb"] for r in rows]
    folds = np.asarray([r["fold"] for r in rows])
    ep = np.asarray([r[pk] for r in rows], float)
    em = np.asarray([r[mk] for r in rows], float)
    base = np.asarray([r[base_key] for r in rows], float)
    bm = em < ep
    orac = np.minimum(ep, em)
    coin = 0.5 * (ep + em)
    fixed = np.where(loo_majority(bm) < 0, em, ep)
    out = {"axis": name, "n": len(rows), "region_mean": float(base.mean()),
           "always_plus": float(ep.mean()), "always_minus": float(em.mean()),
           "coin": float(coin.mean()), "fixed_LOO": float(fixed.mean()),
           "oracle_sign": float(orac.mean()), "frac_minus_better": float(bm.mean()),
           "oracle_vs_fixed": L.report(orac, fixed, folds, pdbs), "arms": {}}
    rng = np.random.default_rng(seed)
    for key in ("z_rg_disto", "z_rg_pool", "z_rg_emit", "z_rg_nat_ORACLE"):
        x = np.asarray([r[key] for r in rows], float)
        sg = loo_sign_rule(x, bm)
        ev = np.where(sg < 0, em, ep)
        acc = float(((sg < 0) == bm).mean())
        # permutation null: shuffle the proxy across targets, refit the LOO rule
        pa, pe = [], []
        for _ in range(n_perm):
            xs = rng.permutation(x)
            s2 = loo_sign_rule(xs, bm)
            pa.append(((s2 < 0) == bm).mean())
            pe.append(np.where(s2 < 0, em, ep).mean())
        out["arms"][key] = {
            "emitted": float(ev.mean()), "sign_accuracy": acc,
            "perm_null_acc_mean": float(np.mean(pa)),
            "perm_p_acc": float((np.asarray(pa) >= acc).mean()),
            "perm_null_emit_mean": float(np.mean(pe)),
            "perm_p_emit": float((np.asarray(pe) <= ev.mean()).mean()),
            "vs_fixed": L.report(ev, fixed, folds, pdbs),
            "vs_oracle": L.report(ev, orac, folds, pdbs),
            "frac_of_oracle_channel": float((fixed.mean() - ev.mean()) /
                                            (fixed.mean() - orac.mean()))
            if fixed.mean() != orac.mean() else float("nan")}
    return out


def main():
    rows = build()
    res = {"n": len(rows), "axes": {}}
    spec = [("p75_plus", "p75_minus", "FULL POOL, top-75", "pool_mean"),
            ("p24_plus", "p24_minus", "FULL POOL, top-24", "pool_mean"),
            ("b75_plus", "b75_minus", "IN-BAND (pool_best+1.5), top-75", "band_mean"),
            ("b24_plus", "b24_minus", "IN-BAND (pool_best+1.5), top-24", "band_mean")]
    for pk, mk, nm, bk in spec:
        res["axes"][nm] = axis(rows, pk, mk, nm, bk)
    res["per_target"] = rows
    L.jwrite("info_signpool", res)

    print("PART B.3 -- DOES THE IN-BAND SIGN RESULT TRANSFER TO THE 126 REAL POOLS?\n")
    for nm, a in res["axes"].items():
        print(f"--- {nm} (n={a['n']}) ---")
        print(f"  region mean {a['region_mean']:.3f}  always+ {a['always_plus']:.3f}  "
              f"always- {a['always_minus']:.3f}  coin {a['coin']:.3f}")
        print(f"  FIXED sign (LOO) {a['fixed_LOO']:.3f}   ORACLE sign {a['oracle_sign']:.3f}"
              f"   frac '-' better {a['frac_minus_better']:.3f}")
        print(f"  ORACLE - fixed : {L.fmt_paired(a['oracle_vs_fixed'])}")
        print(f"  {'proxy':<20}{'emit':>8}{'acc':>7}{'permP':>7}{'ofORC':>7}  vs fixed")
        for k, v in a["arms"].items():
            print(f"  {k:<20}{v['emitted']:>8.3f}{v['sign_accuracy']:>7.3f}"
                  f"{v['perm_p_acc']:>7.3f}{v['frac_of_oracle_channel']:>7.2f}  "
                  f"{L.fmt_paired(v['vs_fixed'])}")
        print()
    return res


if __name__ == "__main__":
    main()
