"""FAIL18 forensics, step 13: SS / short-range-shape as a RETRIEVAL KEY.

E9 showed the fibril subclass is lost at the QUERY: 371 windows per target lie within
2.0 A of the native but BLOSUM ranks the best one at position ~12,000, so only ~17 reach
K=500.  S6-6 tested native SS as a post-hoc FILTER (worth 0.21 A); SS as a retrieval key
-- changing which windows enter the pool at all -- has never been tested.

Keys (each re-ranks the WHOLE universe, then takes K=500):
  blosum        the shipped key: BLOSUM62 sum vs the target          (control)
  shortD        z(blosum) + w * z(-mean |d_cand(i,i+k) - d_pred(i,i+k)|, k=2,3,4)
                d_pred from the SHIPPED distogram -> fully deployable
  ssprop        z(blosum) + w * z(agreement with a propensity-derived SS string)
  ssdisto       z(blosum) + w * z(agreement with an SS string read off the distogram's
                own i,i+4 expected distance)                          -> deployable
  ORACLE_ss     z(blosum) + w * z(agreement with the NATIVE SS string)  [ceiling]
  ORACLE_rg     z(blosum) + w * z(-|rg_cand - rg_native|)               [ceiling]
Metrics: pool best, number of band members in the pool, and top-75 best after the
UNCHANGED shipped distogram filter is applied to the new pool.
"""
from __future__ import annotations
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.environ.setdefault("OMP_NUM_THREADS", "2")
from s12 import instrument as I
from s12.fail_esmrescore import z
from s12.fail_recognise import HELIX_PROP, SHEET_PROP

WGRID = [0.5, 1.0, 2.0]


def ss_from_torsions(PHI, PSI):
    """Cheap Ramachandran-region H/E/C assignment, vectorised over a whole universe."""
    phi = np.degrees(PHI); psi = np.degrees(PSI)
    H = (phi > -160) & (phi < -20) & (psi > -120) & (psi < 50)
    E = (phi > -180) & (phi < -40) & ((psi > 90) | (psi < -150))
    out = np.zeros(phi.shape, np.int8)          # 0 = C
    out[E] = 2
    out[H & ~E] = 1
    return out


def main():
    tg = I.targets()
    rows = []
    for k, t in enumerate(tg):
        pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
        u = I.load_univ(pdb)
        W = u["W"]; rr = u["rr"]; sim = np.asarray(u["sim"], float)
        nat = u["nat_ca"]
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n); sep = j - i
        exp = np.asarray(dg["expected"], float)
        Dp = np.zeros((n, n)); Dp[i, j] = exp; Dp[j, i] = exp

        # short-range shape key
        parts = []
        for kk in (2, 3, 4):
            if n <= kk:
                continue
            a = np.arange(n - kk)
            dc = np.linalg.norm(W[:, a + kk] - W[:, a], axis=-1)
            parts.append(np.abs(dc - Dp[a, a + kk][None, :]).mean(1))
        shortD = -np.mean(parts, 0)

        # SS strings
        SSc = ss_from_torsions(u["PHI"], u["PSI"])              # (nw, n) codes
        # propensity guess
        hp = np.array([HELIX_PROP.get(c, 0.7) for c in seq])
        bp = np.array([SHEET_PROP.get(c, 1.0) for c in seq])
        kern = np.ones(5) / 5.0
        want_prop = np.where(np.convolve(1 - hp / 3.16, kern, "same") >
                             np.convolve(bp / 1.7, kern, "same"), 1, 2).astype(np.int8)
        # distogram-derived: helix if predicted d(i,i+4) < 8 A, strand if > 11 A
        want_dg = np.full(n, 0, np.int8)
        a4 = np.arange(max(n - 4, 0))
        d4 = Dp[a4, a4 + 4] if n > 4 else np.array([])
        for p_ in range(n):
            lo = max(0, p_ - 4); hi = min(len(d4), p_ + 1)
            if hi <= lo:
                continue
            m = d4[lo:hi].mean()
            want_dg[p_] = 1 if m < 8.0 else (2 if m > 11.0 else 0)
        import peptide_db
        pep = peptide_db.by_pdb(pdb)
        nat_ss_str = I.ss_of(pep.phi, pep.psi)
        want_nat = np.array([{"H": 1, "E": 2, "C": 0}[c] for c in nat_ss_str], np.int8)

        def agr(want):
            return (SSc == want[None, :]).mean(1)
        c = W - W.mean(1, keepdims=True)
        rgs = np.sqrt((c ** 2).sum(2).mean(1))
        rgn = float(np.sqrt(((nat - nat.mean(0)) ** 2).sum(1).mean()))

        sigs = dict(shortD=shortD, ssprop=agr(want_prop), ssdisto=agr(want_dg),
                    ORACLE_ss=agr(want_nat), ORACLE_rg=-np.abs(rgs - rgn))
        band_all = rr <= rr.min() + I.BAND
        rec = dict(pdb=pdb, n=n, fold=fold, fail18=pdb in I.FAIL18, arms={})

        def evaluate(p):
            rrp = rr[p]
            Dpair = I.pair_dists(W[p], i, j).astype(np.float32).astype(float)
            sc = I.shipped_score(dg, Dpair)
            sub = np.argsort(sc, kind="stable")[:I.M]
            b = rrp <= rrp.min() + I.BAND
            return dict(pool_best=float(rrp.min()), n_band=int(band_all[p].sum()),
                        top75_best=float(rrp[sub].min()),
                        recall=int(b[sub].any()),
                        top75_mean=float(rrp[sub].mean()))

        rec["arms"]["blosum"] = evaluate(I.pool_idx(u))
        zs = z(sim)
        for name, v in sigs.items():
            for w in WGRID:
                p = np.argsort(-(zs + w * z(v)), kind="stable")[:I.K]
                rec["arms"][f"{name}@{w}"] = evaluate(p)
        rows.append(rec)
        if (k + 1) % 25 == 0:
            print(f"  {k+1}/126 free={I.free_gb():.1f}", flush=True)
    I.write("fail_retrieve", rows)

    f = [r for r in rows if r["fail18"]]; o = [r for r in rows if not r["fail18"]]
    hdr = json.load(open(os.path.join(ROOT, "s12", "results", "fail_headers.json")))["per_target"]
    fib = [r for r in rows if hdr[r["pdb"]]["fibril"]]
    folds = np.array([r["fold"] for r in rows])
    base_pb = np.array([r["arms"]["blosum"]["pool_best"] for r in rows])
    base_t = np.array([r["arms"]["blosum"]["top75_best"] for r in rows])
    print(f"\n== retrieval keys: pool best / top-75 best (ORACLE eval) ==")
    print(f"{'key':16s}{'poolb126':>10s}{'d':>7s}{'ci95':>18s}{'poolbF':>8s}{'poolbFib':>9s}"
          f"{'t75_126':>9s}{'d':>7s}{'t75F':>7s}{'recF':>6s}{'drop10':>8s}")
    summ = {}
    for a in list(rows[0]["arms"]):
        PB = np.array([r["arms"][a]["pool_best"] for r in rows])
        T = np.array([r["arms"][a]["top75_best"] for r in rows])
        st = I.paired(PB, base_pb, folds=folds); stt = I.paired(T, base_t, folds=folds)
        summ[a] = dict(pool_best=float(PB.mean()), top75_best=float(T.mean()),
                       pool_best_fail=float(np.mean([r["arms"][a]["pool_best"] for r in f])),
                       pool_best_fibril=float(np.mean([r["arms"][a]["pool_best"] for r in fib])),
                       top75_best_fail=float(np.mean([r["arms"][a]["top75_best"] for r in f])),
                       recall_fail=float(np.mean([r["arms"][a]["recall"] for r in f])),
                       paired_pool=st, paired_top75=stt)
        s = summ[a]
        print(f"{a:16s}{s['pool_best']:10.3f}{st['mean_diff']:7.3f} "
              f"[{st['ci95'][0]:6.3f},{st['ci95'][1]:6.3f}]{s['pool_best_fail']:8.3f}"
              f"{s['pool_best_fibril']:9.3f}{s['top75_best']:9.3f}{stt['mean_diff']:7.3f}"
              f"{s['top75_best_fail']:7.3f}{s['recall_fail']:6.2f}"
              f"{(stt['drop_top10_mean_diff'] or 0):8.3f}")
    I.write("fail_retrieve_summary", summ)


if __name__ == "__main__":
    main()
