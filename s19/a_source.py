"""SPRINT 19, AGENT A -- WHERE DOES THE COHERENCE COME FROM?  The fork that gates everything.

`a_coh` (n=126) showed that destroying the residual's cross-pair SIGN correlation -- every pair
keeping its own magnitude and its own weight -- moves the deployed refinement 3.610 -> 2.368 A,
129% of the whole Sprint-18 gap, while a PERFECTLY coherent error of matched magnitude is the
worst arm on the board.  The coherence IS the mechanism.  This module asks whether any
predictor change can reach it.

  If the coherent component is SHARED between the deployed MLP, a second seed, a different
  ARCHITECTURE on the same inputs, and the retrieval pool's own mean-distance error, it is an
  inputs / identifiability limit and no loss, architecture or ensemble touches it.
  If it is FAMILY-SPECIFIC, a predictor intervention is worth real money.

THE DECOMPOSITION, and why it is the right one.  Two predictors can agree strongly on the part
of their error that no structure could realise, and that would tell us nothing.  So the error
of each family is split EXACTLY, using that family's own terminal fit X_F:

    r_F        = dhat_F - dtrue                 the raw error
    r_coh_F    = d(X_F) - dtrue                 the part a real structure DOES realise
    r_inc_F    = dhat_F - d(X_F)                the part no structure can realise
    r_F        = r_coh_F + r_inc_F              EXACT, an identity, not a finding

`r_coh_F` is a genuine structure's distance error, so cross-family alignment of `r_coh` is
exactly the question the coordinator asked.

THE NULL FOR "SHARED".  Two independently seeded members of the SAME architecture set the
ceiling for sharing-because-same-model.  Cross-architecture sharing is judged against THAT,
never against zero.  `d20` and `d20_s1` are that pair (identical architecture and training,
`torch.manual_seed` 0 vs 1).

THE MOST DIRECT STATISTIC.  The coherent component is a structure.  So ask whether the families
converge on the SAME WRONG STRUCTURE: CA-RMSD(X_A, X_B) against each family's own distance to
the native.  Every family is fitted from the IDENTICAL start, so agreement could be shared
start inertia -- the mandatory control is CA-RMSD(X_A, start) and CA-RMSD(X_A, X_A^signflip),
the structure the SAME family reaches once its own sign coherence is destroyed.

SURGERY (through the same deployed fit; ORACLE; magnitude-matched twins).
  minus_shared / only_shared   the deployed residual with its component along PairNet's
                               residual regressed out, and only that component
  minus_pool   / only_pool     the same against the retrieval pool's mean-distance error

Run:  python -m s19.a_source     (requires s19/cache/zoo_probs.npz, built by a_models)
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s19 import a_lib as L                   # noqa: E402
from s19 import a_fit as F                   # noqa: E402
from s19 import a_models as M                # noqa: E402
from s12 import instrument as I              # noqa: E402
from s14 import avgspace as AV               # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import seed as SD                   # noqa: E402

OUT = os.path.join(L.RESULTS, "a_source.json")
FAMS = ["deployed", "d20", "d20_s1", "pairnet", "poolmean", "sepprior", "helix"]
PAIRS = [("d20", "d20_s1", "SAME architecture, different seed -- the CEILING for sharing"),
         ("deployed", "d20", "same architecture, different regularisation"),
         ("deployed", "pairnet", "DIFFERENT architecture, same inputs"),
         ("d20", "pairnet", "DIFFERENT architecture, matched regularisation"),
         ("deployed", "poolmean", "not a network at all -- retrieval"),
         ("pairnet", "poolmean", "cross-architecture vs retrieval"),
         ("deployed", "sepprior", "ZERO-INFORMATION: sequence-blind separation prior"),
         ("deployed", "helix", "ZERO-INFORMATION: a constant ideal alpha-helix"),
         ("pairnet", "sepprior", "zero-information vs the joint architecture")]
ARMS = ["real", "signflip",
        "minus_shared", "only_shared", "minus_pool", "only_pool",
        "minus_sep", "only_sep",
        "minus_shared_m", "only_shared_m", "minus_pool_m", "only_pool_m",
        "minus_sep_m", "only_sep_m"]


def _corr(a, b, m=None):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if m is not None:
        a, b = a[m], b[m]
    if len(a) < 4:
        return np.nan
    a = a - a.mean()
    b = b - b.mean()
    den = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / den) if den > 0 else np.nan


def _sign_coh(r, i, j):
    """sign-product over pair-pairs SHARING a residue minus over pairs sharing none."""
    s = np.sign(np.asarray(r, float))
    P = len(s)
    if P < 5:
        return np.nan
    share = (i[:, None] == i[None, :]) | (i[:, None] == j[None, :]) | \
            (j[:, None] == i[None, :]) | (j[:, None] == j[None, :])
    tri = np.triu(np.ones((P, P), bool), 1)
    pr = s[:, None] * s[None, :]
    a = pr[share & tri]
    b = pr[(~share) & tri]
    if not len(a) or not len(b):
        return np.nan
    return float(a.mean() - b.mean())


def _fit_struct(dhat, sd, i, j, phi0, psi0):
    t = np.maximum(np.asarray(dhat, float), 2.0)
    ph, ps, f = A.fit(t, sd, i, j, phi0, psi0)
    ca = I.build_ca(ph, ps)
    return ca, np.linalg.norm(ca[i] - ca[j], axis=1)


def run(tg=None):
    tg = tg if tg is not None else I.targets()
    probs = M.probs_for(tg)
    fams = M.build_data(tg, probs)
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)

    fam = {}
    for t in tg:
        p = t["pdb"]
        base = dict(fams["deployed"][p])
        W = np.asarray(AV.top75_windows(p)[0], float)
        D = np.linalg.norm(W[:, base["i"], :] - W[:, base["j"], :], axis=-1)
        base["dhat"] = D.mean(0)
        base["sd"] = np.maximum(D.std(0), 1e-3)
        fam[p] = base
    deb = L.debias_fns(fam, pdbs, folds)
    for p in pdbs:
        fam[p]["dhat"] = np.maximum(fam[p]["dhat"] - deb[fam[p]["fold"]](fam[p]["sep"]), 2.0)
    fams["poolmean"] = fam

    # ---- ZERO-INFORMATION FAMILY 1: the sequence-blind separation prior.
    # d(sep) = mean TRUE CA-CA distance at that separation over the TRAINING folds only, so
    # it is leave-fold-out and carries no information about the target beyond its length.
    # Its `sd` is the training-fold spread at that separation, so it is a real predictor.
    base0 = fams["deployed"]
    sp_all = {f: {} for f in sorted(set(folds.tolist()))}
    for f in sp_all:
        tr = [p for p in pdbs if base0[p]["fold"] != f]
        S = np.concatenate([base0[p]["sep"] for p in tr])
        T = np.concatenate([base0[p]["dtrue"] for p in tr])
        for s in np.unique(S):
            m = S == s
            sp_all[f][int(s)] = (float(T[m].mean()), float(max(T[m].std(), 0.1)))
    fam = {}
    for p in pdbs:
        b = dict(base0[p])
        tbl = sp_all[b["fold"]]
        mu = np.array([tbl.get(int(s), (10.0, 3.0))[0] for s in b["sep"]])
        sg = np.array([tbl.get(int(s), (10.0, 3.0))[1] for s in b["sep"]])
        b["dhat"] = np.maximum(mu, 2.0)
        b["sd"] = np.maximum(sg, 1e-3)
        fam[p] = b
    fams["sepprior"] = fam

    # ---- ZERO-INFORMATION FAMILY 2: a constant ideal alpha-helix (brief section 8's
    # mandated plausible-but-uninformative reference measure).  Its `sd` is the deployed
    # distogram's, so the fit it produces is comparable.
    fam = {}
    for p in pdbs:
        b = dict(base0[p])
        n = b["n"]
        ca = I.build_ca(np.full(n, np.deg2rad(-63.0)), np.full(n, np.deg2rad(-42.0)))
        b["dhat"] = np.maximum(np.linalg.norm(ca[b["i"]] - ca[b["j"]], axis=1), 2.0)
        fam[p] = b
    fams["helix"] = fam

    data, _, _ = L.gather_all(tg)
    rows = []
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        d = data[pdb]
        i, j, sd, nat, n = d["i"], d["j"], d["sd"], d["nat"], d["n"]
        dtrue = d["dtrue"]
        phi0, psi0, avg = F.start(pdb, d["seq"], d["fold"])
        start_ca = I.build_ca(phi0, psi0)
        rng = SD.stable_rng(pdb, "s19A_source")
        uni = (d["nmode"] <= 1)
        conf = (sd <= np.median(sd))
        uc = uni & conf

        R, RC, RI, X = {}, {}, {}, {}
        for f in FAMS:
            dh = fams[f][pdb]["dhat"]
            # each family is fitted with ITS OWN weights -- that is the object it produces
            ca, dfit = _fit_struct(dh, fams[f][pdb]["sd"], i, j, phi0, psi0)
            X[f] = ca
            R[f] = dh - dtrue
            RC[f] = dfit - dtrue                      # the realisable part -- EXACT identity
            RI[f] = dh - dfit                         # the unrealisable remainder
        r0 = R["deployed"]
        P = len(r0)

        e = {"pdb": pdb, "n": n, "fold": d["fold"], "npairs": P,
             "frac_uc": float(uc.mean()), "resid_rms": float(np.sqrt((r0 ** 2).mean())),
             "start_rmsd": float(I.ca_rmsd(start_ca, nat))}
        for f in FAMS:
            e[f"rmsd|{f}"] = float(I.ca_rmsd(X[f], nat))
            e[f"tostart|{f}"] = float(I.ca_rmsd(X[f], start_ca))
            e[f"coh|{f}"] = _sign_coh(R[f], i, j)
            e[f"cohUC|{f}"] = _sign_coh(R[f][uc], i[uc], j[uc]) if uc.sum() > 5 else np.nan
            e[f"cohC|{f}"] = _sign_coh(RC[f], i, j)
            e[f"share_coh|{f}"] = float((RC[f] ** 2).sum() / max((R[f] ** 2).sum(), 1e-12))
        e["coh|null"] = _sign_coh(rng.permutation(r0), i, j)

        for a, b, _lab in PAIRS:
            e[f"corr|{a}|{b}"] = _corr(R[a], R[b])
            e[f"corrUC|{a}|{b}"] = _corr(R[a], R[b], uc)
            e[f"corrC|{a}|{b}"] = _corr(RC[a], RC[b])
            e[f"corrI|{a}|{b}"] = _corr(RI[a], RI[b])
            e[f"xrmsd|{a}|{b}"] = float(I.ca_rmsd(X[a], X[b]))

        # the same-family whitened structure -- the null for structural agreement
        sf = dtrue + r0 * rng.choice([-1.0, 1.0], size=P)
        ca_sf, _ = _fit_struct(sf, sd, i, j, phi0, psi0)
        e["xrmsd|deployed|signflip"] = float(I.ca_rmsd(X["deployed"], ca_sf))
        e["rmsd|signflip"] = float(I.ca_rmsd(ca_sf, nat))

        # ---------------- surgery
        def split(other):
            x = R[other]
            xc = x - x.mean()
            b = float((xc * (r0 - r0.mean())).sum() / max((xc * xc).sum(), 1e-12))
            comp = b * xc + r0.mean()
            return comp, r0 - comp

        sh_only, sh_minus = split("pairnet")
        pl_only, pl_minus = split("poolmean")
        sp_only, sp_minus = split("sepprior")
        rms = float(np.sqrt((r0 ** 2).mean()))
        fields = {"real": r0, "signflip": r0 * rng.choice([-1.0, 1.0], size=P),
                  "minus_shared": sh_minus, "only_shared": sh_only,
                  "minus_pool": pl_minus, "only_pool": pl_only,
                  "minus_sep": sp_minus, "only_sep": sp_only}
        for k in ("minus_shared", "only_shared", "minus_pool", "only_pool",
                  "minus_sep", "only_sep"):
            v = fields[k]
            fields[k + "_m"] = v * (rms / max(float(np.sqrt((v ** 2).mean())), 1e-9))
        for a in ARMS:
            fld = np.maximum(dtrue + fields[a], 2.0)
            e[a] = F.fit_rmsd(fld, sd, i, j, phi0, psi0, nat)[0]
            e[a + "_rms"] = float(np.sqrt((fields[a] ** 2).mean()))
        rows.append(e)
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"),
                      default=float)
    json.dump({"rows": rows, "complete": len(rows) == len(tg)}, open(OUT, "w"), default=float)
    if len(rows) == len(tg):
        open(os.path.join(L.RESULTS, "a_source.COMPLETE"), "w").write("ok\n")
    report(rows)
    return rows


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    rng = SD.stable_rng("s19A", "source", "report")
    g = lambda k: np.array([r[k] for r in rows], float)               # noqa: E731
    nan = np.nanmean
    tab = {}
    print(f"\n=== SOURCE ATTRIBUTION   n = {len(rows)} ===")
    print(f"  unimodal AND confident pairs: {g('frac_uc').mean():.1%} of all pairs")
    print(f"  shared start (projected coordinate average): {g('start_rmsd').mean():.3f} A\n")

    print("  --- each family's own fit, and how much of its error is REALISABLE ---")
    print(f"  {'family':<12}{'RMSD':>8}{'to start':>10}{'coherent share of SSE':>24}"
          f"{'sign-coherence':>16}")
    for f in FAMS:
        print(f"  {f:<12}{g('rmsd|'+f).mean():>8.3f}{g('tostart|'+f).mean():>10.3f}"
              f"{nan(g('share_coh|'+f)):>24.3f}{nan(g('coh|'+f)):>16.4f}")
    print(f"  {'PERM null':<12}{'':>8}{'':>10}{'':>24}{nan(g('coh|null')):>16.4f}")
    tab["families"] = {f: {"rmsd": float(g("rmsd|" + f).mean()),
                           "to_start": float(g("tostart|" + f).mean()),
                           "coherent_share": float(nan(g("share_coh|" + f))),
                           "sign_coh": float(nan(g("coh|" + f))),
                           "sign_coh_uniconf": float(nan(g("cohUC|" + f)))} for f in FAMS}
    tab["sign_coh_null"] = float(nan(g("coh|null")))

    print("\n  --- sign-coherence on UNIMODAL + CONFIDENT pairs only (the coordinator's ask) ---")
    for f in FAMS:
        print(f"    {f:<12}all {nan(g('coh|'+f)):+.4f}   unimodal+confident "
              f"{nan(g('cohUC|'+f)):+.4f}")

    print("\n  --- CROSS-FAMILY ALIGNMENT OF THE ERROR ---")
    print(f"  {'pair':<24}{'raw r':>9}{'COHERENT':>10}{'incoherent':>12}"
          f"{'uni+conf':>10}{'RMSD(X_A,X_B)':>15}   what it is")
    tab["pairs"] = {}
    for a, b, lab in PAIRS:
        print(f"  {a+'~'+b:<24}{nan(g(f'corr|{a}|{b}')):>9.3f}"
              f"{nan(g(f'corrC|{a}|{b}')):>10.3f}{nan(g(f'corrI|{a}|{b}')):>12.3f}"
              f"{nan(g(f'corrUC|{a}|{b}')):>10.3f}{g(f'xrmsd|{a}|{b}').mean():>15.3f}   {lab}")
        tab["pairs"][f"{a}~{b}"] = {
            "raw": float(nan(g(f"corr|{a}|{b}"))), "coherent": float(nan(g(f"corrC|{a}|{b}"))),
            "incoherent": float(nan(g(f"corrI|{a}|{b}"))),
            "uni_conf": float(nan(g(f"corrUC|{a}|{b}"))),
            "xrmsd": float(g(f"xrmsd|{a}|{b}").mean()), "label": lab}
    print(f"  {'deployed~ITS OWN signflip':<24}{'':>9}{'':>10}{'':>12}{'':>10}"
          f"{g('xrmsd|deployed|signflip').mean():>15.3f}   THE STRUCTURAL NULL")
    tab["structural_null"] = float(g("xrmsd|deployed|signflip").mean())

    print("\n  --- SURGERY THROUGH THE DEPLOYED FIT (ORACLE) ---")
    real = g("real")
    print(f"  {'arm':<18}{'RMSD':>8}{'med':>8}{'residRMS':>10}   vs real")
    tab["arms"] = {}
    for a in ARMS:
        v = g(a)
        d = v - real
        m, lo, hi = L.boot(d, rng)
        print(f"  {a:<18}{v.mean():>8.3f}{np.median(v):>8.3f}{g(a+'_rms').mean():>10.3f}"
              f"   {m:+.3f} [{lo:+.3f},{hi:+.3f}]  {int((d<0).sum())}W/{int((d>0).sum())}L")
        tab["arms"][a] = {"rmsd": float(v.mean()), "diff": m, "ci": [lo, hi],
                          "W": int((d < 0).sum()), "L": int((d > 0).sum())}
    json.dump(tab, open(os.path.join(L.RESULTS, "a_source_report.json"), "w"), indent=1,
              default=float)
    return tab


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
