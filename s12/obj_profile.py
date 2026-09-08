"""The SEPARATION-PROFILE channel: measure it, and try to predict it better (all LFO).

Established by `obj_iso2` / `obj_arm4b`: essentially all of the usable information in the
CA-CA distance objective is the per-|i-j| mean distance profile q_s.  This module asks how
well the shipped distogram predicts that profile, and whether anything deployable predicts
it better.  Arms:

  disto        the shipped distogram's own profile (baseline)
  pool         the K=500 pool's mean profile (native-free typicality, no model at all)
  debias       disto minus the per-shell mean bias estimated on the OTHER folds
  cal          per-shell across-target affine recalibration q* ~ a_s + b_s q_disto  (LFO)
  cal2         as `cal` but with the pool profile as a second regressor            (LFO)
  expand       disto expanded about the training-fold shell mean by a factor k
  ridge        LFO ridge from ESM-2 PCA-32 mean/std pooled features -> profile
  ridge_all    LFO ridge from [ESM pooled, composition, length, disto profile, pool profile]
  oracle       the true profile (DIAGNOSTIC ceiling)

The scored objective is always  target_p = q'_{sep(p)} + (E[d]_p - q_disto_{sep(p)}) , i.e.
only the profile is swapped; the distogram's pair-specific residual is left untouched.
"""
from __future__ import annotations
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from s12 import instrument as I
from s12 import obj_common as OC
from s12 import obj_arm4 as A4

SMAX = 15


def profiles(d):
    """(shell ids, true profile, distogram profile, pool profile) for one target."""
    sep = d["sep"]; us = np.unique(sep)
    qt = np.array([d["dtrue"][sep == s].mean() for s in us])
    qp = np.array([d["exp"][sep == s].mean() for s in us])
    qo = np.array([d["D"][:, sep == s].mean() for s in us])
    return us, qt, qp, qo


def esm_feat(seq):
    from s12 import esm_bank
    B = esm_bank.load()
    if seq not in B:
        return None
    e = B[seq][0]                                   # (n, 32)
    return np.concatenate([e.mean(0), e.std(0), [len(seq)]])


def comp_feat(seq):
    A = I.ALPHABET
    c = np.array([seq.count(a) for a in A], float) / len(seq)
    from core.predict import _PROP_TABLE
    idx = np.array([A.index(a) for a in seq])
    P = _PROP_TABLE[idx]
    return np.concatenate([c, P.mean(0), P.std(0), [len(seq)]])


def build():
    tg = OC.targets()
    rec = []
    for t in tg:
        d = OC.load(t["pdb"])
        us, qt, qp, qo = profiles(d)
        f_esm = esm_feat(t["seq"])
        rec.append(dict(pdb=t["pdb"], fold=t["fold"], n=t["n"], seq=t["seq"], us=us,
                        qt=qt, qp=qp, qo=qo, comp=comp_feat(t["seq"]),
                        esm=f_esm if f_esm is not None else np.zeros(65)))
    return rec


def ridge_fit(X, Y, lam=1.0):
    X = np.concatenate([X, np.ones((len(X), 1))], 1)
    A = X.T @ X + lam * np.eye(X.shape[1]); A[-1, -1] -= lam
    return np.linalg.solve(A, X.T @ Y)


def ridge_pred(W, X):
    return np.concatenate([X, np.ones((len(X), 1))], 1) @ W


def fit_all(rec, lam=3.0):
    """Per-shell LFO fits.  Returns per target a dict of predicted profiles by arm."""
    folds = sorted({r["fold"] for r in rec})
    out = {r["pdb"]: {} for r in rec}
    for f in folds:
        tr = [r for r in rec if r["fold"] != f]
        te = [r for r in rec if r["fold"] == f]
        # ---- per-shell scalar fits
        bias, cal, cal2 = {}, {}, {}
        for s in range(2, SMAX + 1):
            xs = np.array([r["qp"][list(r["us"]).index(s)] for r in tr if s in r["us"]])
            ys = np.array([r["qt"][list(r["us"]).index(s)] for r in tr if s in r["us"]])
            zs = np.array([r["qo"][list(r["us"]).index(s)] for r in tr if s in r["us"]])
            if len(xs) < 12:
                bias[s] = 0.0; cal[s] = (0.0, 1.0); cal2[s] = None; continue
            bias[s] = float((xs - ys).mean())
            b = np.polyfit(xs, ys, 1); cal[s] = (float(b[1]), float(b[0]))
            M = np.stack([xs, zs, np.ones_like(xs)], 1)
            cal2[s] = np.linalg.lstsq(M, ys, rcond=None)[0]
        # ---- vector fits (fixed-length profile s=2..8, pad with the shell mean)
        SS = list(range(2, 9))
        def vec(r, key):
            u = list(r["us"])
            return np.array([r[key][u.index(s)] if s in u else np.nan for s in SS])
        Ytr = np.array([vec(r, "qt") for r in tr]); Xd = np.array([vec(r, "qp") for r in tr])
        Xo = np.array([vec(r, "qo") for r in tr])
        mu = np.nanmean(Ytr, 0)
        fill = lambda M: np.where(np.isnan(M), mu[None, :], M)
        Ytr, Xd, Xo = fill(Ytr), fill(Xd), fill(Xo)
        Etr = np.array([r["esm"] for r in tr]); Ctr = np.array([r["comp"] for r in tr])
        Fsets = {"ridge": np.concatenate([Etr, Ctr], 1),
                 "ridge_all": np.concatenate([Etr, Ctr, Xd, Xo], 1)}
        Ws = {}
        for k, F in Fsets.items():
            m, sd = F.mean(0), F.std(0) + 1e-6
            Ws[k] = (m, sd, ridge_fit((F - m) / sd, Ytr, lam))
        for r in te:
            u = list(r["us"]); qp = r["qp"]; qo = r["qo"]
            o = {}
            o["disto"] = qp.copy()
            o["pool"] = qo.copy()
            o["oracle"] = r["qt"].copy()
            o["debias"] = np.array([qp[k] - bias.get(s, 0.0) for k, s in enumerate(u)])
            o["cal"] = np.array([cal.get(s, (0.0, 1.0))[0] + cal.get(s, (0.0, 1.0))[1] * qp[k]
                                 for k, s in enumerate(u)])
            o["cal2"] = np.array([(cal2[s] @ np.array([qp[k], qo[k], 1.0])) if cal2.get(s) is not None
                                  else qp[k] for k, s in enumerate(u)])
            for kk in ("ridge", "ridge_all"):
                m, sd, W = Ws[kk]
                F = np.concatenate([r["esm"], r["comp"]] +
                                   ([vec(r, "qp"), vec(r, "qo")] if kk == "ridge_all" else []))
                F = np.where(np.isnan(F), 0.0, F)
                p = ridge_pred(W, ((F - m) / sd)[None, :])[0]
                o[kk] = np.array([p[SS.index(s)] if s in SS else qp[k] for k, s in enumerate(u)])
            for kexp in (1.3, 1.7, 2.2):
                mu_s = {s: float(np.mean([rr["qp"][list(rr["us"]).index(s)] for rr in tr if s in rr["us"]]))
                        for s in u if any(s in rr["us"] for rr in tr)}
                o[f"expand{kexp}"] = np.array([mu_s.get(s, qp[k]) + kexp * (qp[k] - mu_s.get(s, qp[k]))
                                               for k, s in enumerate(u)])
            out[r["pdb"]] = o
    return out


def evaluate(prof, names=None):
    tg = OC.targets()
    names = names or sorted(next(iter(prof.values())).keys())
    rows = []
    for t in tg:
        d = OC.load(t["pdb"]); sep = d["sep"]; us = np.unique(sep)
        qp = np.array([d["exp"][sep == s].mean() for s in us])
        base = d["exp"] - np.array([qp[list(us).index(s)] for s in sep])   # pair-specific residual
        qt = np.array([d["dtrue"][sep == s].mean() for s in us])
        for nm in names:
            q = prof[t["pdb"]][nm]
            tgt = base + np.array([q[list(us).index(s)] for s in sep])
            sc = OC.score_l1(d["D"], tgt)
            e = OC.emit(d, sc, lam=None)
            rows.append(dict(pdb=t["pdb"], fold=t["fold"], tag=nm, avg_rmsd=e["avg_rmsd"],
                             sub_best=e["sub_best"], argmin_rmsd=e["argmin_rmsd"],
                             prof_mae=float(np.abs(q - qt).mean()),
                             mae=float(np.abs(tgt - d["dtrue"]).mean()),
                             rho=A4.spearman(sc, d["rr"])))
    return rows


def report(rows):
    import collections
    by = collections.defaultdict(list)
    for r in rows:
        by[r["tag"]].append(r)
    hdr = f"{'profile arm':14s} {'profMAE':>8s} {'pairMAE':>8s} {'avg':>7s} {'FAIL18':>7s} {'oth108':>7s} {'top75b':>7s} {'rho':>6s}"
    print(hdr); print("-" * len(hdr))
    for tag, v in sorted(by.items(), key=lambda kv: np.mean([x["avg_rmsd"] for x in kv[1]])):
        f = lambda q, sel=None: float(np.nanmean([x[q] for x in v if sel is None or sel(x)]))
        print(f"{tag:14s} {f('prof_mae'):8.3f} {f('mae'):8.3f} {f('avg_rmsd'):7.3f} "
              f"{f('avg_rmsd', lambda x: x['pdb'] in I.FAIL18):7.3f} "
              f"{f('avg_rmsd', lambda x: x['pdb'] not in I.FAIL18):7.3f} {f('sub_best'):7.3f} {f('rho'):6.3f}")


if __name__ == "__main__":
    rec = build()
    prof = fit_all(rec)
    rows = evaluate(prof)
    I.write("obj_profile", {"rows": rows})
    report(rows)
