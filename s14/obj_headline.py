"""SPRINT 14, OBJ -- the headline leave-fold-out run: arms, ablations, nulls, curve.

Every arm is trained on targets whose PINNED FOLD differs from the test target's, with
per-target centering and standardised ridge, and evaluated on a fixed uniform
262,144-configuration evaluation set drawn from the test target's full enumeration.

Arms
  learned    the model: distance-binned pair potential + Rama locals + globals +
             composition x global interactions, terminal residues masked
  nomask     identical but WITHOUT masking the two CA-RMSD-inert terminal residues
  noint      identical but WITHOUT the composition x global interaction block
  band       trained on the near-native 0.1% RMSD band instead of uniform samples
  seqblind   NULL: no residue identity anywhere
  permuted   NULL: labels permuted within each training target
  randfeat   NULL: features replaced by random sparse noise of matched density
  leak       POSITIVE CONTROL: the true within-target RMSD percentile as an extra input
             column at train AND test time.  Not a model; a loudness calibrator.

Baselines: Legacy, the 1-local prior, and random (0 by construction on the delta axes).

    python -m s14.obj_headline arms
    python -m s14.obj_headline curve
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
from s14 import obj_floor as F                     # noqa: E402

LAM = 0.03            # PREDEFINED from the first exploratory sweep; fixed for every arm
NUNI = 8000
INTER = slice(M.NPAIRF + M.NLOCF + M.NGLOB, M.DIM - 1)

FZ = {}


def fz_for(mask_terminal=True):
    if mask_terminal not in FZ:
        FZ[mask_terminal] = M.Featurizer(cache_dir=True, mask_terminal=mask_terminal)
    return FZ[mask_terminal]


ARMSPEC = {
    "learned":  dict(),
    "nomask":   dict(mask_terminal=False),
    "noint":    dict(drop=INTER),
    "band":     dict(n_uniform=0, n_band=NUNI, band_frac=0.001),
    "seqblind": dict(sequence_blind=True),
    "permuted": dict(permute_label=True),
    "randw":    dict(),          # post-fit: weights randomised, see fit_arm
    "leak":     dict(leak=True),
}


class FoldBank:
    """Per-fold normal-equation blocks for one feature kind, built in ONE pass.

    Normal equations are additive across per-target centered blocks, so the
    leave-fold-out training block is `TOTAL - fold`.  That turns 5 folds x 4 feature
    kinds = 20 featurization passes into 4, which is the difference between this run
    finishing and not.

    Three arms are derived algebraically from the same pass rather than refitted:
      noint     drop the interaction columns  ->  zero those rows/cols of XtX and Xty
      permuted  labels permuted within target ->  same XtX, a second Xty accumulated
      leak      append the label as a column  ->  the augmented system is
                [[XtX, Xty], [Xty', y'y]] with rhs [Xty, y'y], all already known
    """

    def __init__(self, dim=None):
        self.dim = dim or M.DIM
        self.XtX = np.zeros((self.dim, self.dim))
        self.Xty = np.zeros(self.dim)
        self.Xtyp = np.zeros(self.dim)
        self.yty = 0.0
        self.n = 0

    def add_target(self, X, y, yperm):
        N = X.shape[0]
        mu = np.asarray(X.sum(0)).ravel() / N
        ybar, ypbar = float(y.mean()), float(yperm.mean())
        self.XtX += (X.T @ X).toarray() - N * np.outer(mu, mu)
        self.Xty += X.T @ y - N * mu * ybar
        self.Xtyp += X.T @ yperm - N * mu * ypbar
        self.yty += float(y @ y) - N * ybar * ybar
        self.n += N

    def __isub__(self, o):
        self.XtX -= o.XtX; self.Xty -= o.Xty; self.Xtyp -= o.Xtyp
        self.yty -= o.yty; self.n -= o.n
        return self

    def copy(self):
        o = FoldBank(self.dim)
        o.XtX = self.XtX.copy(); o.Xty = self.Xty.copy(); o.Xtyp = self.Xtyp.copy()
        o.yty = self.yty; o.n = self.n
        return o

    def _solve(self, A, b, lam):
        N = max(self.n, 1)
        A = A / N; b = b / N
        s = np.sqrt(np.maximum(np.diag(A), 0.0)); s[s < 1e-9] = 1.0
        D = 1.0 / s
        A = A * D[:, None] * D[None, :] + lam * np.eye(len(A))
        # Do NOT exempt the intercept from the penalty here: per-target centering makes
        # the intercept column exactly zero, so exempting it leaves a zero row AND column
        # and the system is singular.  It only ever shifts every score by a constant, so
        # rankings are unaffected either way -- but leaving it penalised keeps the solve
        # conditioned and the weight norms interpretable.
        return np.linalg.solve(A, b * D) * D

    def w_learned(self, lam=None):
        return self._solve(self.XtX.copy(), self.Xty.copy(), lam or LAM)

    def w_permuted(self, lam=None):
        return self._solve(self.XtX.copy(), self.Xtyp.copy(), lam or LAM)

    def w_noint(self, lam=None):
        A = self.XtX.copy(); b = self.Xty.copy()
        A[INTER, :] = 0.0; A[:, INTER] = 0.0; b[INTER] = 0.0
        return self._solve(A, b, lam or LAM)

    def w_leak(self, lam=None):
        d = self.dim
        A = np.zeros((d + 1, d + 1)); bb = np.zeros(d + 1)
        A[:d, :d] = self.XtX
        A[:d, d] = A[d, :d] = self.Xty
        A[d, d] = self.yty
        bb[:d] = self.Xty; bb[d] = self.yty
        o = FoldBank(d + 1); o.n = self.n
        return o._solve(A, bb, lam or LAM)


def build_bank(targets, folds, kind, seed=0):
    """One featurization pass -> (TOTAL, {fold: FoldBank})."""
    fz = fz_for(kind != "nomask")
    rng = np.random.default_rng(seed)
    prng = np.random.default_rng(seed + 17)
    per = {}
    total = FoldBank()
    for p in targets:
        en = T.Enum(p)
        idx = (T.sample_idx(en, 0, NUNI, rng, band_frac=0.001) if kind == "band"
               else T.sample_idx(en, NUNI, 0, rng))
        X = fz.csr(p, T._decode(idx, en.n, en.k), sequence_blind=(kind == "seqblind"))
        y = en.pct[idx].astype(np.float64)
        f = folds[p]
        if f not in per:
            per[f] = FoldBank()
        yp = prng.permutation(y)      # the SAME permutation must enter both, or the
        per[f].add_target(X, y, yp)   # leave-fold-out subtraction mixes the two labels
        total.add_target(X, y, yp)
        del X
    return total, per


_XCACHE = {}


def blocks_for(train, kind, seed=0):
    """CSR design blocks for one training set, built ONCE and shared by every arm that
    uses the same features and the same sample.  Five of the nine arms differ only in
    the label or in a post-fit weight transform, so this is a ~3x saving over refitting
    each arm from raw features."""
    key = (tuple(sorted(train)), kind, seed)
    if key in _XCACHE:
        return _XCACHE[key]
    rng = np.random.default_rng(seed)
    fz = fz_for(kind != "nomask")
    out = []
    for p in train:
        en = T.Enum(p)
        if kind == "band":
            idx = T.sample_idx(en, 0, NUNI, rng, band_frac=0.001)
        else:
            idx = T.sample_idx(en, NUNI, 0, rng)
        S = T._decode(idx, en.n, en.k)
        X = fz.csr(p, S, sequence_blind=(kind == "seqblind"))
        out.append((p, idx, X, en.pct[idx].astype(np.float64)))
    _XCACHE.clear()                 # one training set live at a time; these are large
    _XCACHE[key] = out
    return out


def fit_blocks(blocks, drop=None, permute=False, leak=False, seed=0):
    from scipy import sparse
    acc = M.RidgeAccumulator(M.DIM + (1 if leak else 0))
    rng = np.random.default_rng(seed + 17)
    for _p, _idx, X, y in blocks:
        if permute:
            y = rng.permutation(y)
        if drop is not None:
            keep = np.ones(X.shape[1]); keep[drop] = 0.0
            X = X @ sparse.diags(keep)
        if leak:
            X = sparse.hstack([X, sparse.csr_matrix(y.reshape(-1, 1))], format="csr")
        acc.add(X, y, center=True)
    return acc.solve(LAM, standardize=True)


def fit_all(train, arms, seed=0):
    """Fit every requested arm on `train`, sharing design matrices. -> {arm: (w, mt)}"""
    out = {}
    need = set(arms)
    shared = need & {"learned", "noint", "permuted", "randw", "leak"}
    if shared:
        B = blocks_for(train, "uniform", seed)
        if "learned" in need or "randw" in need:
            out["learned"] = (fit_blocks(B, seed=seed), True)
        if "noint" in need:
            out["noint"] = (fit_blocks(B, drop=INTER, seed=seed), True)
        if "permuted" in need:
            out["permuted"] = (fit_blocks(B, permute=True, seed=seed), True)
        if "leak" in need:
            out["leak"] = (fit_blocks(B, leak=True, seed=seed), True)
        if "randw" in need:
            w = out["learned"][0]
            rng = np.random.default_rng(seed + 4242)
            rw = np.zeros_like(w)
            for a, b in [(0, M.NPAIRF), (M.NPAIRF, M.NPAIRF + M.NLOCF),
                         (M.NPAIRF + M.NLOCF, M.NPAIRF + M.NLOCF + M.NGLOB),
                         (M.NPAIRF + M.NLOCF + M.NGLOB, M.DIM - 1)]:
                rw[a:b] = rng.standard_normal(b - a) * w[a:b].std()
            out["randw"] = (rw, True)
        if "learned" not in arms:
            out.pop("learned", None)
    for a, kind, mt in (("seqblind", "seqblind", True), ("nomask", "nomask", False),
                        ("band", "band", True)):
        if a in need:
            out[a] = (fit_blocks(blocks_for(train, kind, seed), seed=seed), mt)
    return out


def fit_arm(train, arm, seed=0):
    if arm == "randw":
        # NULL: a random weight vector whose per-feature-family scale matches the fitted
        # model's.  This is the sharp null for THIS feature space, because an early
        # diagnostic showed the highest-variance columns dominate any weight vector --
        # so "some linear functional of these features orders structures" is not
        # evidence that the FITTED one learned anything.
        w, fz, mt = fit_arm(train, "learned", seed=seed)
        rng = np.random.default_rng(seed + 4242)
        out = np.zeros_like(w)
        blocks = [(0, M.NPAIRF), (M.NPAIRF, M.NPAIRF + M.NLOCF),
                  (M.NPAIRF + M.NLOCF, M.NPAIRF + M.NLOCF + M.NGLOB),
                  (M.NPAIRF + M.NLOCF + M.NGLOB, M.DIM - 1)]
        for a, b in blocks:
            out[a:b] = rng.standard_normal(b - a) * w[a:b].std()
        return out, fz, mt
    s = dict(ARMSPEC[arm])
    mt = s.pop("mask_terminal", True)
    bf = s.pop("band_frac", None)
    fz = fz_for(mt)
    kw = dict(n_uniform=s.pop("n_uniform", NUNI), n_band=s.pop("n_band", 0))
    if bf is not None:
        nb = kw["n_band"]
        kw = dict(n_uniform=0, n_band=0)
        s["idx_fn"] = lambda en, rng: T.sample_idx(en, 0, nb, rng, band_frac=bf)
    w = T.fit(train, fz, lam=LAM, center=True, seed=seed, **kw, **s)
    return w, fz, mt


def score_arm(w, arm, pdb_id, idx, mt=True):
    fz = fz_for(mt)
    en = T.Enum(pdb_id)
    S = T._decode(idx, en.n, en.k)
    ww = np.asarray(w)[:M.DIM].copy()
    if ARMSPEC[arm].get("drop") is not None:
        ww[ARMSPEC[arm]["drop"]] = 0.0
    e = fz.score_with(ww, pdb_id, S, ARMSPEC[arm].get("sequence_blind", False))
    if ARMSPEC[arm].get("leak"):
        e = e + float(np.asarray(w)[M.DIM]) * en.pct[idx].astype(np.float64)
    return e


def stats_of(e, r, seed=0):
    o = np.argsort(e, kind="mergesort")
    dec = o[:max(50, len(e) // 10)]
    rng = np.random.default_rng(seed)
    pi = F.pair_discrimination(e, r, rng, npair=300_000, sub=dec)
    pg = F.pair_discrimination(e, r, rng, npair=300_000)
    lo = e.min()
    tie = np.flatnonzero(e == lo)                 # tie-averaged argmin (project finding)
    return dict(rho_global=spearman(e, r), rho_decile=spearman(e[dec], r[dec]),
                d_decile=float(r[dec].mean() - r.mean()),
                d_top100=float(r[o[:100]].mean() - r.mean()),
                d_top1000=float(r[o[:1000]].mean() - r.mean()),
                d_argmin=float(r[tie].mean() - r.mean()),
                argmin_ties=int(len(tie)),
                inband_lt15=float(np.nanmean([pi[k][0] for k in
                                              ("0.25-0.5", "0.5-1.0", "1.0-1.5")
                                              if k in pi])),
                global_lt15=float(np.nanmean([pg[k][0] for k in
                                              ("0.25-0.5", "0.5-1.0", "1.0-1.5")
                                              if k in pg])),
                space_mean=float(r.mean()))


FIELDS = ["rho_global", "rho_decile", "inband_lt15", "global_lt15",
          "d_decile", "d_top100", "d_argmin"]


KINDS = {"uniform": ["learned", "noint", "permuted", "leak", "randw"],
         "seqblind": ["seqblind"], "nomask": ["nomask"], "band": ["band"]}


def run_banked(targets, arms=None, seed=0):
    arms = arms or list(ARMSPEC)
    folds = {p: T.Enum(p).fold for p in targets}
    fl = sorted(set(folds.values()))
    print(f"{len(targets)} targets, folds "
          f"{ {f: sum(1 for p in targets if folds[p]==f) for f in fl} }", flush=True)
    evx = {p: X.eval_idx(T.Enum(p)) for p in targets}
    rr = {p: T.Enum(p).rmsd[evx[p]].astype(np.float64) for p in targets}
    res = {a: {} for a in arms}
    base = {"legacy": {}, "prior": {}}
    t0 = time.time()
    for p in targets:
        en = T.Enum(p)
        base["legacy"][p] = stats_of(en.legacy[evx[p]].astype(np.float64), rr[p])
        base["prior"][p] = stats_of(en.prior[evx[p]].astype(np.float64), rr[p])
    print(f"baselines [{time.time()-t0:.0f}s]", flush=True)

    for kind, kin_arms in KINDS.items():
        want = [a for a in kin_arms if a in arms]
        if not want:
            continue
        total, per = build_bank(targets, folds, kind, seed=seed)
        print(f"bank[{kind}] built [{time.time()-t0:.0f}s]", flush=True)
        for f in fl:
            tr = total.copy()
            tr -= per[f]
            test = [p for p in targets if folds[p] == f]
            ws = {}
            if "learned" in want or "randw" in want:
                ws["learned"] = (tr.w_learned(), True)
            for a, fn, mt in (("noint", tr.w_noint, True),
                              ("permuted", tr.w_permuted, True),
                              ("leak", tr.w_leak, True),
                              ("seqblind", tr.w_learned, True),
                              ("nomask", tr.w_learned, False),
                              ("band", tr.w_learned, True)):
                if a in want:
                    ws[a] = (fn(), mt)
            if "randw" in want:
                w = ws["learned"][0]
                rng = np.random.default_rng(seed + 4242)
                rw = np.zeros_like(w)
                for a0, b0 in [(0, M.NPAIRF), (M.NPAIRF, M.NPAIRF + M.NLOCF),
                               (M.NPAIRF + M.NLOCF, M.NPAIRF + M.NLOCF + M.NGLOB),
                               (M.NPAIRF + M.NLOCF + M.NGLOB, M.DIM - 1)]:
                    rw[a0:b0] = rng.standard_normal(b0 - a0) * w[a0:b0].std()
                ws["randw"] = (rw, True)
            for a in want:
                w, mt = ws[a]
                for p in test:
                    res[a][p] = stats_of(score_arm(w, a, p, evx[p], mt), rr[p])
                print(f"  fold {f} {a:9s} "
                      + " ".join(f"{p}:{res[a][p]['d_decile']:+.2f}" for p in test)
                      + f"  [{time.time()-t0:.0f}s]", flush=True)
            del tr, ws
        del total, per
    summarise(res, base, targets, folds)
    I.write("s14_obj_headline",
            {"lam": LAM, "n_uniform": NUNI, "targets": targets,
             "folds": {p: int(v) for p, v in folds.items()},
             "arms": res, "baselines": base}, n_expected=len(targets))
    return res, base


def run(targets, arms=None, seed=0):
    arms = arms or list(ARMSPEC)
    folds = {p: T.Enum(p).fold for p in targets}
    print(f"{len(targets)} targets, folds "
          f"{ {f: sum(1 for p in targets if folds[p]==f) for f in sorted(set(folds.values()))} }",
          flush=True)
    evx = {p: X.eval_idx(T.Enum(p)) for p in targets}
    rr = {p: T.Enum(p).rmsd[evx[p]].astype(np.float64) for p in targets}

    res = {a: {} for a in arms}
    base = {"legacy": {}, "prior": {}}
    t0 = time.time()
    for p in targets:
        en = T.Enum(p)
        base["legacy"][p] = stats_of(en.legacy[evx[p]].astype(np.float64), rr[p])
        base["prior"][p] = stats_of(en.prior[evx[p]].astype(np.float64), rr[p])
    print(f"baselines [{time.time()-t0:.0f}s]", flush=True)

    for f in sorted(set(folds.values())):
        test = [p for p in targets if folds[p] == f]
        train = [p for p in targets if folds[p] != f]
        raise RuntimeError("superseded by run_banked")
    summarise(res, base, targets, folds)
    I.write("s14_obj_headline",
            {"lam": LAM, "n_uniform": NUNI, "targets": targets,
             "folds": {p: int(v) for p, v in folds.items()},
             "arms": {a: {p: v for p, v in d.items()} for a, d in res.items()},
             "baselines": base}, n_expected=len(targets))
    return res, base


def summarise(res, base, targets, folds):
    fa = np.array([folds[p] for p in targets])
    print(f"\n{'arm':10s} " + " ".join(f"{f:>12s}" for f in FIELDS) + "   W/L(dec)")
    allr = list(base.items()) + list(res.items())
    for nm, d in allr:
        v = {f: np.array([d[p][f] for p in targets], float) for f in FIELDS}
        print(f"{nm:10s} " + " ".join(f"{np.nanmean(v[f]):12.3f}" for f in FIELDS)
              + f"   {int((v['d_decile']<0).sum())}/{int((v['d_decile']>0).sum())}")
    print("\n-- paired vs RANDOM (0) on d_decile; negative = better than random --")
    for nm, d in allr:
        a = np.array([d[p]["d_decile"] for p in targets])
        pr = I.paired(a, np.zeros_like(a), folds=fa, names=targets)
        print(f"{nm:10s} {pr['mean_diff']:+.3f} CI[{pr['ci95'][0]:+.3f},"
              f"{pr['ci95'][1]:+.3f}] W/L={pr['n_better']}/{pr['n_worse']} "
              f"drop10={pr['drop_top10_mean_diff']:+.3f} "
              f"folds={ {k: round(v,3) for k,v in pr['per_fold'].items()} }")
    print("\n-- paired vs RANDOM (0) on d_top100 --")
    for nm, d in allr:
        a = np.array([d[p]["d_top100"] for p in targets])
        pr = I.paired(a, np.zeros_like(a), folds=fa, names=targets)
        print(f"{nm:10s} {pr['mean_diff']:+.3f} CI[{pr['ci95'][0]:+.3f},"
              f"{pr['ci95'][1]:+.3f}] W/L={pr['n_better']}/{pr['n_worse']} "
              f"drop10={pr['drop_top10_mean_diff']:+.3f}")
    print("\n-- paired: each arm vs LEGACY on d_decile (negative = arm better) --")
    lg = np.array([base["legacy"][p]["d_decile"] for p in targets])
    for nm, d in res.items():
        a = np.array([d[p]["d_decile"] for p in targets])
        pr = I.paired(a, lg, folds=fa, names=targets)
        print(f"{nm:10s} {pr['mean_diff']:+.3f} CI[{pr['ci95'][0]:+.3f},"
              f"{pr['ci95'][1]:+.3f}] W/L={pr['n_better']}/{pr['n_worse']}")


def curve(targets, sizes=(1, 2, 4, 8, 16), nrep=3, seed=0):
    """Learning curve in number of training targets, with the leaked-label control.

    Test set is one whole pinned fold; the training pool is every target in the other
    folds, so no size ever contains a fold-mate of a test target.
    """
    folds = {p: T.Enum(p).fold for p in targets}
    tf = max(set(folds.values()), key=lambda f: sum(1 for p in targets if folds[p] == f))
    test = [p for p in targets if folds[p] == tf]
    pool = [p for p in targets if folds[p] != tf]
    print(f"held-out fold {tf}: {test}\npool ({len(pool)}): {pool}", flush=True)
    evx = {p: X.eval_idx(T.Enum(p)) for p in test}
    rr = {p: T.Enum(p).rmsd[evx[p]].astype(np.float64) for p in test}
    sizes = [m for m in sizes if m <= len(pool)]
    if len(pool) not in sizes:
        sizes.append(len(pool))
    out = {}
    t0 = time.time()
    for arm in ("learned", "leak"):
        for m in sizes:
            acc = []
            for rep in range(nrep if m < len(pool) else 1):
                rng = np.random.default_rng(1000 * m + rep + seed)
                tr = list(rng.choice(pool, m, replace=False))
                w, mt = fit_all(tr, [arm], seed=seed)[arm]
                for p in test:
                    acc.append(stats_of(score_arm(w, arm, p, evx[p], mt), rr[p]))
            out[f"{arm}|{m}"] = {f: float(np.nanmean([a[f] for a in acc]))
                                 for f in FIELDS}
            out[f"{arm}|{m}"]["n_obs"] = len(acc)
            print(f"{arm:8s} m={m:2d} " + "  ".join(
                f"{f}={out[f'{arm}|{m}'][f]:+.3f}" for f in FIELDS)
                + f"  [{time.time()-t0:.0f}s]", flush=True)
    I.write("s14_obj_curve", {"test_fold": int(tf), "test": test, "pool": pool,
                              "rows": out}, n_expected=1)
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "arms"
    tg = [a for a in sys.argv[2:] if not a.startswith("-")]
    tg = tg or [t["pdb"] for t in I.targets() if E.have(t["pdb"])]
    if mode == "arms":
        run_banked(tg)
    elif mode == "curve":
        curve(tg)
