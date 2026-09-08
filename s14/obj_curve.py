"""SPRINT 14, OBJ -- the learning curve, its leaked-label positive control, and the
out-of-distribution (generator-shift) test.

Sprint 12's decisive negative was a FLAT learning curve against a leaked-label control
that was already loud at n = 8.  This runs the identical diagnostic on the torsion-space
objective.  If the curve is flat while the leak control climbs, the limit is SIGNAL, not
sample size, and no amount of extra enumerated targets or model capacity will help.

The leak arm appends the TRUE within-target RMSD percentile as one extra input column at
BOTH train and test time.  It is not a model, it is a loudness calibrator: it says what a
working curve looks like in this harness.

Generator shift: the model is fit on configurations from one generator and evaluated on
configurations from another (uniform-random vs sampled from the leakage-safe 1-local
prior).  A VQE generates its own distribution, so an objective that only ranks the
distribution it was fitted on is useless for VQE.

    python -m s14.obj_curve curve
    python -m s14.obj_curve ood
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                    # noqa: E402
from s13.qarch_lib import spearman                 # noqa: E402
from s14 import obj_enum as E                      # noqa: E402
from s14 import obj_model as M                     # noqa: E402
from s14 import obj_train as T                     # noqa: E402
from s14 import obj_exp as X                       # noqa: E402

LAM = 0.03
SIZES = [1, 2, 4, 8]
NREP = 3


def score_arm(fz, w, pdb_id, idx, leak=False, sequence_blind=False):
    en = T.Enum(pdb_id)
    S = T._decode(idx, en.n, en.k)
    e = fz.score_with(np.asarray(w)[:M.DIM], pdb_id, S, sequence_blind)
    if leak:
        e = e + float(np.asarray(w)[M.DIM]) * en.pct[idx].astype(np.float64)
    return e


def stats(e, r, seed=0):
    o = np.argsort(e, kind="mergesort")
    dec = o[:max(50, len(e) // 10)]
    from s14 import obj_floor as F
    pi = F.pair_discrimination(e, r, np.random.default_rng(seed), npair=200_000, sub=dec)
    return dict(rho_global=spearman(e, r), rho_decile=spearman(e[dec], r[dec]),
                d_decile=float(r[dec].mean() - r.mean()),
                d_top100=float(r[o[:100]].mean() - r.mean()),
                inband_lt15=float(np.nanmean(
                    [pi[k][0] for k in ("0.25-0.5", "0.5-1.0", "1.0-1.5") if k in pi])))


def curve(targets, test=None, sizes=SIZES, nrep=NREP, seed=0, n_uniform=20000):
    fz = M.Featurizer(cache_dir=True)
    rng0 = np.random.default_rng(seed)
    test = test or list(rng0.choice(targets, 4, replace=False))
    pool = [p for p in targets if p not in test]
    print(f"test targets: {test}\npool: {pool} (max m={len(pool)})", flush=True)
    sizes = [m for m in sizes if m <= len(pool)] + ([len(pool)] if len(pool) not in sizes
                                                    else [])
    evx = {p: X.eval_idx(T.Enum(p)) for p in test}
    rr = {p: T.Enum(p).rmsd[evx[p]].astype(np.float64) for p in test}
    out = {}
    t0 = time.time()
    for arm, leak in (("learned", False), ("leak_control", True)):
        for m in sizes:
            acc = []
            for rep in range(nrep if m < len(pool) else 1):
                rng = np.random.default_rng(1000 * m + rep + seed)
                tr = list(rng.choice(pool, m, replace=False))
                w = T.fit(tr, fz, n_uniform=n_uniform, lam=LAM, center=True, leak=leak,
                          seed=seed)
                for p in test:
                    e = score_arm(fz, w, p, evx[p], leak=leak)
                    acc.append(stats(e, rr[p]))
            k = (arm, m)
            out[k] = {f: float(np.nanmean([a[f] for a in acc])) for f in acc[0]}
            out[k]["n_obs"] = len(acc)
            print(f"{arm:14s} m={m:2d}  " + "  ".join(
                f"{f}={out[k][f]:+.3f}" for f in ("rho_global", "rho_decile", "d_decile",
                                                  "d_top100", "inband_lt15"))
                  + f"   [{time.time()-t0:.0f}s]", flush=True)
    I.write("s14_obj_curve", {"test": test, "pool": pool, "sizes": sizes,
                              "rows": {f"{a}|{m}": v for (a, m), v in out.items()}},
            n_expected=1)
    return out


def ood(targets, seed=0, n_uniform=20000, n_eval=120_000):
    """Fit on one generator's configurations, evaluate on another's."""
    fz = M.Featurizer(cache_dir=True)
    rng = np.random.default_rng(seed)

    def unif_idx(en, r):
        return r.choice(en.B, min(n_uniform, en.B), replace=False)

    def prior_idx(en, r):
        return np.unique(T.sample_prior_idx(en, n_uniform, r))

    gens = {"uniform": unif_idx, "prior": prior_idx}
    print(f"{'fit_on':9s} {'eval_on':9s} {'rho_g':>8s} {'rho_dec':>8s} "
          f"{'d_decile':>9s} {'d_top100':>9s} {'inband<1.5':>11s}")
    res = {}
    for gtr, ftr in gens.items():
        rows = {g: [] for g in gens}
        for p in targets:
            tr = [q for q in targets if q != p]
            w = T.fit(tr, fz, lam=LAM, center=True, seed=seed,
                      idx_fn=lambda en, r: ftr(en, r))
            en = T.Enum(p)
            for gte, fte in gens.items():
                idx = np.unique(fte(en, np.random.default_rng(seed + 7)))
                if len(idx) > n_eval:
                    idx = idx[:n_eval]
                e = score_arm(fz, w, p, idx)
                rows[gte].append(stats(e, en.rmsd[idx].astype(np.float64)))
        for gte, acc in rows.items():
            v = {f: float(np.nanmean([a[f] for a in acc])) for f in acc[0]}
            res[f"{gtr}->{gte}"] = v
            print(f"{gtr:9s} {gte:9s} {v['rho_global']:+8.3f} {v['rho_decile']:+8.3f} "
                  f"{v['d_decile']:+9.3f} {v['d_top100']:+9.3f} {v['inband_lt15']:11.3f}",
                  flush=True)
    I.write("s14_obj_ood", {"rows": res}, n_expected=1)
    return res


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "curve"
    tg = [a for a in sys.argv[2:] if not a.startswith("-")]
    tg = tg or [t["pdb"] for t in I.targets() if E.have(t["pdb"])]
    if mode == "curve":
        curve(tg)
    elif mode == "ood":
        ood(tg)
