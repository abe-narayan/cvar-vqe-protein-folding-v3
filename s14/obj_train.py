"""SPRINT 14, OBJ, STEPS 3-4 -- train the learned objective and try to break it.

Training label: within-target RMSD percentile (rank / (B-1)) on the fully enumerated
k=4 state space.  Ridge on streamed normal equations.  Leave-fold-out on the pinned
folds; the pilot mode is leave-one-target-out and is labelled as such.

Evaluation is always on the FULL enumeration of a held-out target, and always reports
objective quality and structural quality separately.

    python -m s14.obj_train pilot            # leave-one-target-out on what is enumerated
    python -m s14.obj_train lfo              # leave-fold-out, the headline
    python -m s14.obj_train curve            # learning curve + leaked-label control
    python -m s14.obj_train ood              # generator-shift robustness
    python -m s14.obj_train nulls            # sequence-blind / permuted / random-feature
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                       # noqa: E402
from s13.qarch_lib import spearman                    # noqa: E402
from s14 import obj_enum as E                         # noqa: E402
from s14 import obj_model as M                        # noqa: E402
from s14 import obj_floor as F                        # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
os.makedirs(RESULTS, exist_ok=True)

LAMS = [3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0]


# ------------------------------------------------------------------ data access
class Enum:
    """Cached view of one target's enumeration: rmsd, its percentile, the baselines."""

    _cache = {}

    def __new__(cls, pdb_id):
        if pdb_id in cls._cache:
            return cls._cache[pdb_id]
        self = super().__new__(cls)
        z = E.load_enum(pdb_id)
        self.pdb = pdb_id
        self.n = int(z["n"])
        self.k = int(z["k"])
        self.seq = str(z["seq"])
        self.fold = int(z["fold"])
        self.rmsd = np.asarray(z["rmsd"], np.float32)
        self.legacy = np.asarray(z["legacy"], np.float32)
        self.prior = np.asarray(z["prior"], np.float32)
        self.snap = int(z["snap_index"])
        self.B = len(self.rmsd)
        o = np.argsort(self.rmsd, kind="mergesort")
        pct = np.empty(self.B, np.float32)
        pct[o] = np.arange(self.B, dtype=np.float32) / (self.B - 1)
        self.pct = pct
        cls._cache[pdb_id] = self
        return self

    def states(self, idx):
        return E.enum_configs(self.n, self.k, 0, 1) if idx is None else \
            _decode(np.asarray(idx, np.int64), self.n, self.k)


def _decode(r, n, k):
    out = np.empty((len(r), n), np.int8)
    for i in range(n):
        out[:, i] = (r // (k ** (n - 1 - i))) % k
    return out


def sample_idx(en, n_uniform, n_band, rng, band_frac=0.01):
    parts = []
    if n_uniform:
        parts.append(rng.choice(en.B, min(n_uniform, en.B), replace=False))
    if n_band:
        nb = max(2000, int(en.B * band_frac))
        band = np.argpartition(en.rmsd, nb)[:nb]
        parts.append(rng.choice(band, min(n_band, len(band)), replace=False))
    return np.unique(np.concatenate(parts))


def sample_prior_idx(en, m, rng):
    """Configurations drawn from the leakage-safe 1-local prior (a different generator)."""
    from s14 import obj_model as MM
    fz = MM.Featurizer(k=en.k, cache_dir=True)
    P = fz.prior_table(en.pdb)
    cum = np.cumsum(P, axis=1)
    u = rng.random((m, en.n))
    S = (u[:, :, None] > cum[None, :, :-1]).sum(2)
    return np.ravel_multi_index(tuple(S.T), (en.k,) * en.n)


# ------------------------------------------------------------------ fit / eval
def fit(train_pdbs, fz, n_uniform=30000, n_band=0, seed=0, lam=None,
        sequence_blind=False, permute_label=False, random_feature=False,
        leak=False, idx_fn=None, verbose=False, center=True, label="pct",
        drop=None, fz_override=None):
    """Accumulate normal equations over `train_pdbs` and return {lam: w} or w.

    `center=True` fits within-target deviations (see RidgeAccumulator.add).
    `label`: "pct" = within-target RMSD percentile; "rmsd" = raw A.
    """
    acc = M.RidgeAccumulator(M.DIM + (1 if leak else 0))
    rng = np.random.default_rng(seed)
    for p in train_pdbs:
        en = Enum(p)
        idx = (idx_fn(en, rng) if idx_fn is not None
               else sample_idx(en, n_uniform, n_band, rng))
        S = _decode(idx, en.n, en.k)
        X = fz.csr(p, S, sequence_blind=sequence_blind)
        y = (en.pct[idx] if label == "pct" else en.rmsd[idx]).astype(np.float64)
        if permute_label:
            y = rng.permutation(y)
        if random_feature:
            from scipy import sparse
            Xr = sparse.random(X.shape[0], X.shape[1], density=min(1.0, X.nnz / X.size),
                               random_state=int(rng.integers(1 << 30)),
                               format="csr", dtype=np.float32)
            X = Xr
        if drop is not None:
            from scipy import sparse
            keep = np.ones(X.shape[1]); keep[drop] = 0.0
            X = X @ sparse.diags(keep)
        if leak:
            from scipy import sparse
            X = sparse.hstack([X, sparse.csr_matrix(y.reshape(-1, 1))], format="csr")
        acc.add(X, y, center=center)
        if verbose:
            print(f"    fit {p}: {len(idx)} rows", flush=True)
    if lam is not None:
        return acc.solve(lam)
    return acc.solve_path(LAMS)


def eval_target(pdb_id, w, fz, sequence_blind=False, rng=None, full=True,
                sub_idx=None, extra_dim=0):
    """All Step-2 statistics for the learned objective on one target."""
    en = Enum(pdb_id)
    rng = rng or np.random.default_rng(0)
    ww = np.asarray(w)[:M.DIM] if extra_dim else np.asarray(w)
    if full and sub_idx is None:
        S = _decode(np.arange(en.B, dtype=np.int64), en.n, en.k)
        e = fz.score_with(ww, pdb_id, S, sequence_blind)
        r = en.rmsd.astype(np.float64)
        snap = en.snap
    else:
        idx = sub_idx
        S = _decode(idx, en.n, en.k)
        e = fz.score_with(ww, pdb_id, S, sequence_blind)
        r = en.rmsd[idx].astype(np.float64)
        snap = int(np.where(idx == en.snap)[0][0]) if (idx == en.snap).any() else None
    st = F.analyse_column(e, r, snap, rng, "learned")
    st["pdb"] = pdb_id
    st["fold"] = en.fold
    st["rmsd_mean_all"] = float(en.rmsd.mean())
    return st


def baselines(pdb_id, rng=None):
    en = Enum(pdb_id)
    rng = rng or np.random.default_rng(0)
    r = en.rmsd.astype(np.float64)
    out = {}
    for nm, e in (("legacy", en.legacy.astype(np.float64)),
                  ("prior", en.prior.astype(np.float64))):
        out[nm] = F.analyse_column(e, r, en.snap, rng, nm)
    out["random"] = dict(rmsd_argmin=float(r.mean()), rmsd_decile_mean=float(r.mean()),
                         top100_mean_rmsd=float(r.mean()), rho_global=0.0,
                         rho_decile=0.0, top1_mean_rmsd=float(r.mean()))
    out["oracle_best"] = dict(rmsd_argmin=float(r.min()))
    return out


HEAD = ["rho_global", "rho_decile", "rmsd_argmin", "top100_mean_rmsd",
        "rmsd_decile_mean", "snap_pct"]


def report(rows, tag, extra=()):
    print(f"\n== {tag} ==")
    print(f"{'pdb':6s} {'fold':>4s} " + " ".join(f"{h:>17s}" for h in HEAD))
    for s in rows:
        print(f"{s['pdb']:6s} {s['fold']:4d} " +
              " ".join(f"{s.get(h, float('nan')):17.3f}" for h in HEAD))
    print(f"{'MEAN':6s} {'':4s} " +
          " ".join(f"{np.nanmean([s.get(h, np.nan) for s in rows]):17.3f}"
                   for h in HEAD))


def inband_pair(stats, upto=1.5):
    """Mean in-band pairwise-ordering accuracy for separations below `upto` A."""
    keys = [f"{lo}-{hi}" for lo, hi in F.DBINS if hi <= upto]
    v = [stats["pair_inband"][kk][0] for kk in keys
         if kk in stats["pair_inband"] and np.isfinite(stats["pair_inband"][kk][0])]
    return float(np.mean(v)) if v else float("nan")


# ------------------------------------------------------------------ experiments
def enumerated():
    return [t["pdb"] for t in I.targets() if E.have(t["pdb"])]


def pilot(pdbs=None, n_uniform=30000, n_band=0, lam=None):
    pdbs = pdbs or enumerated()
    fz = M.Featurizer(cache_dir=True)
    rows, base = [], []
    t0 = time.time()
    for p in pdbs:
        tr = [q for q in pdbs if q != p]
        ws = fit(tr, fz, n_uniform=n_uniform, n_band=n_band)
        # inner selection: pick lambda on the training targets themselves (no test leak)
        if lam is None:
            best, bl = -9, None
            for l, w in ws.items():
                sc = np.mean([spearman(fz.score_with(
                    w, q, _decode(sample_idx(Enum(q), 4000, 0,
                                             np.random.default_rng(7)), Enum(q).n, 4)),
                    Enum(q).rmsd[sample_idx(Enum(q), 4000, 0,
                                            np.random.default_rng(7))])
                    for q in tr[:4]])
                if sc > best:
                    best, bl = sc, l
        else:
            bl = lam
        w = ws[bl]
        s = eval_target(p, w, fz)
        s["lam"] = bl
        rows.append(s)
        base.append(baselines(p))
        print(f"{p}: learned rho_g={s['rho_global']:+.3f} rho_dec={s['rho_decile']:+.3f} "
              f"dec_mean={s['rmsd_decile_mean']:.3f} argmin={s['rmsd_argmin']:.3f} "
              f"(legacy dec {base[-1]['legacy']['rmsd_decile_mean']:.3f}, "
              f"random {base[-1]['random']['rmsd_argmin']:.3f}) lam={bl} "
              f"[{time.time()-t0:.0f}s]", flush=True)
    return rows, base


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "pilot"
    if mode == "pilot":
        rows, base = pilot()
        report(rows, "LEARNED (leave-one-target-out pilot)")
        report([b["legacy"] | {"pdb": r["pdb"], "fold": r["fold"]}
                for r, b in zip(rows, base)], "LEGACY")
        I.write("s14_obj_pilot", {"rows": rows}, n_expected=len(rows))
