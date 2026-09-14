#!/usr/bin/env python
"""s27/run_pool.py -- compute every channel on every target's shipped pool, cache it, and
evaluate every registered configuration on the POINT-CLOUD endpoint (`s27/PREREG.md` H1 to H6).

Per target this writes `s27/cache/<pdb>.npz` (every channel as a (500,) vector plus DIS, LEG,
AMB from the S25 caches) and appends its rows to `s27/results/pool_rows.jsonl`; the run is
resumable.  ORACLE reads happen only inside `oracle_rmsd_of_set` / `oracle_diag`.

Usage:  python s27/run_pool.py [--limit N]
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import d_harness as H             # noqa: E402
from s25 import phys_lib as P              # noqa: E402
from s27 import ham_lib as HL              # noqa: E402

CACHE = os.path.join(HERE, "cache")
RESULTS = os.path.join(HERE, "results")
ROWS = os.path.join(RESULTS, "pool_rows.jsonl")
M = 75
RAMA = np.load(os.path.join(ROOT, "s8", "generate_rama.npz"))["cnt"]

NEW = list(HL.CHANNELS)
LEGT = ["LEG_" + t for t in HL.LEG_TERMS]
SINGLES = NEW + LEGT + ["DIS", "LEG", "AMB"]
COMPLEMENTS = NEW + LEGT + ["LEG", "AMB"]          # each paired with DIS
COMPOSITES = {
    "ROSETTA_LIKE": ["ENV", "CONTACT", "EXVOL", "RG_LAW"],
    "PHYSICS_COMP": ["DSSPHB", "ELEC", "HP", "EXVOL"],
    "CONSIST_COMP": ["CONS", "DMAP_CONS", "TORS_CONS", "POOLGO"],
    "STAT_COMP": ["CONTACT", "DISTPOT", "ENV", "CAGEO"],
}
REJECT_Q = 0.10


def rng_for(pdb, tag):
    """A stable per-(target, purpose) RNG: sha256, not Python's per-process string hash."""
    import hashlib
    h = hashlib.sha256(f"s27|{pdb}|{tag}".encode()).digest()
    return np.random.default_rng(int.from_bytes(h[:8], "little"))


def topm(E, m, key):
    """Top-m of E with EXACT ties broken by the random key, never by array order."""
    E = np.asarray(E, float)
    order = np.lexsort((key, E))
    return order[:m]


def zr(x):
    return HL.zrank(x)


def combine(ch, names, weights=None):
    s = 0.0
    for k, nm in enumerate(names):
        w = 1.0 if weights is None else float(weights[k])
        s = s + w * zr(ch[nm])
    return zr(s)


def oracle_rmsd_of_set(cand, idx):
    return P.rmsd_of_set(cand, idx)


def oracle_diag(cand, E, dis_top):
    """ORACLE diagnostics of one energy vector: rho with the candidate RMSD, in-band rho,
    top-75 overlap with DIS."""
    from scipy.stats import spearmanr
    rr = np.asarray(cand.oracle_rr, float)
    E = np.asarray(E, float).copy()
    if not np.isfinite(E).all():                      # staged rejects carry +inf: rank them last
        fin = np.isfinite(E)
        E[~fin] = (E[fin].max() if fin.any() else 0.0) + 1.0
    rho = float(spearmanr(E, rr).correlation) if np.std(E) > 0 else 0.0
    inb = rr < 3.0
    rho_in = float(spearmanr(E[inb], rr[inb]).correlation) if inb.sum() > 5 and np.std(E[inb]) > 0 else float("nan")
    return dict(rho_pool=rho, rho_inband=rho_in, n_inband=int(inb.sum()))


def channels_for(pdb):
    f = os.path.join(CACHE, f"{pdb}.npz")
    cand, ref = P.channels(pdb)                       # DIS / LEG / AMB with the cache asserts
    if os.path.exists(f):
        z = np.load(f)
        ch = {k: np.asarray(z[k], float) for k in z.files if k != "cost_ms"}
        cost = json.loads(str(z["cost_ms"]))
        return cand, ch, cost
    u = I.load_univ(pdb)
    dg = I.distogram(pdb, cand.seq, cand.fold)
    cx = HL.Context(cand, u, dg, RAMA[cand.fold])
    ch, cost = {}, {}
    for name, fn in HL.CHANNELS.items():
        t0 = time.perf_counter()
        ch[name] = np.asarray(fn(cx), float)
        cost[name] = (time.perf_counter() - t0) * 1e3 / cand.k
    t0 = time.perf_counter()
    lt = HL.legacy_terms(cx)
    dt = (time.perf_counter() - t0) * 1e3 / cand.k
    for k, v in lt.items():
        ch[k] = v
        cost[k] = dt
    ch.update({"DIS": ref["DIS"], "LEG": ref["LEG"], "AMB": ref["AMB"]})
    os.makedirs(CACHE, exist_ok=True)
    np.savez_compressed(f + ".tmp.npz", cost_ms=json.dumps(cost), **ch)
    os.replace(f + ".tmp.npz", f)
    return cand, ch, cost


def evaluate_target(pdb):
    cand, ch, cost = channels_for(pdb)
    k = cand.k
    key = rng_for(pdb, "tiekey").random(k)
    perm_rng = rng_for(pdb, "perm")
    perm = {c: perm_rng.permutation(k) for c in ch}
    rows = []
    dis_top = topm(ch["DIS"], M, key)
    r_dis = oracle_rmsd_of_set(cand, dis_top)
    # random-75 null, 16 draws
    rr = rng_for(pdb, "rand75")
    rnd = np.array([oracle_rmsd_of_set(cand, rr.choice(k, size=M, replace=False)) for _ in range(16)])
    base = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), rmsd_dis75=r_dis,
                rand_mean=float(rnd.mean()), rand_sd=float(rnd.std(ddof=1)),
                oracle_pool_best=float(np.min(cand.oracle_rr)), oracle_pool_mean=float(np.mean(cand.oracle_rr)))

    def add(config, kind, E, extra=None):
        E = np.asarray(E, float)
        top = topm(E, M, key)
        r = oracle_rmsd_of_set(cand, top)
        d = oracle_diag(cand, E, dis_top)
        row = dict(base, config=config, kind=kind, rmsd=r, overlap_dis=float(len(set(top) & set(dis_top)) / M),
                   n_distinct=int(len(np.unique(E))), tie_frac_top=float(1 - len(np.unique(E[top])) / M),
                   **d)
        if extra:
            row.update(extra)
        rows.append(row)

    # H1 singles (+ rho with DIS, cost)
    from scipy.stats import spearmanr
    for c in SINGLES:
        rho_dis = float(spearmanr(ch[c], ch["DIS"]).correlation) if np.std(ch[c]) > 0 else 0.0
        add(c, "single", ch[c], dict(rho_dis=rho_dis, cost_ms=float(cost.get(c, float("nan")))))
    # H2 equal-weight complements, with the permuted control
    for c in COMPLEMENTS:
        add(f"DIS+{c}", "pair", combine(ch, ["DIS", c]))
        chp = dict(ch); chp[c] = ch[c][perm[c]]
        add(f"DIS+{c}~perm", "pair_perm", combine(chp, ["DIS", c]))
    # H3 regulariser weights
    for c in COMPLEMENTS:
        for w in (0.25, 0.5):
            add(f"DIS+{w}*{c}", "weighted", combine(ch, ["DIS", c], [1.0, w]))
    # H4 staged reject
    q = int(round(REJECT_Q * k))
    for c in COMPLEMENTS:
        keep = topm(ch[c], k - q, key)                       # drop the worst q by c
        Ek = np.full(k, np.inf); Ek[keep] = ch["DIS"][keep]
        add(f"REJ[{c}]->DIS", "staged", Ek)
    keep = rng_for(pdb, "rejrand").choice(k, size=k - q, replace=False)
    Ek = np.full(k, np.inf); Ek[keep] = ch["DIS"][keep]
    add("REJ[random]->DIS", "staged_ctrl", Ek)
    # H5 composites, alone and with DIS
    for nm, parts in COMPOSITES.items():
        add(nm, "composite", combine(ch, parts))
        add(f"DIS+{nm}", "pair", combine(ch, ["DIS"] + parts, [len(parts)] + [1.0] * len(parts)))
    # H6 adaptive mixture: weight from the distogram's mean per-pair entropy (native-free)
    dg = I.distogram(pdb, cand.seq, cand.fold)
    pr = np.clip(np.asarray(dg["prob"], float), 1e-12, 1)
    Hbar = float(-(pr * np.log(pr)).sum(1).mean())            # nats, per pair
    hmax = math.log(17.0)
    w_t = 0.5 * Hbar / hmax
    for c in ("CONS", "DSSPHB", "CONTACT", "DISTPOT", "RAMA", "LEG", "AMB"):
        add(f"DIS+adapt({c})", "adaptive", combine(ch, ["DIS", c], [1.0, w_t]), dict(w_t=w_t, entropy_nats=Hbar))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    done = set()
    if os.path.exists(ROWS):
        with open(ROWS, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["pdb"])
                except Exception:
                    pass
    pdbs = P.targets()
    if a.limit:
        pdbs = pdbs[:a.limit]
    t0 = time.time()
    for n, pdb in enumerate(pdbs):
        if pdb in done:
            continue
        t1 = time.time()
        rows = evaluate_target(pdb)
        with open(ROWS, "a", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        print(f"  [{n+1}/{len(pdbs)}] {pdb} n={rows[0]['n']} rows={len(rows)} "
              f"dis75={rows[0]['rmsd_dis75']:.3f} rand={rows[0]['rand_mean']:.3f} {time.time()-t1:.1f}s "
              f"(elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", ROWS)


if __name__ == "__main__":
    main()
