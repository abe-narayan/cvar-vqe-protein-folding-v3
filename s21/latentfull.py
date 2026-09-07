"""s21/latentfull.py -- ENUMERATE THE ENTIRE VQE LATENT SPACE AND TAKE THE EXACT ARGMIN.

WHY THIS CLOSES SOMETHING SIX INSTRUMENTS ONLY BOUNDED.  From source (`s19/qb_lib.py:198`),
`draw_from_basins` takes bits of shape (B, n) with n = RESIDUES and bits in {0,1} -- **one qubit
per residue, two basins each** -- so the entire latent the VQE explores is `2**n`.  Measured over
the 126 tuning targets: min 512, **median 8192**, max 65,536.  At Sprint 20's 8192-evaluation
budget the budget **equals or exceeds the whole latent on 75 of 126 targets**.

    So on most targets the deployed VQE was not searching a space it could not enumerate.
    It was RESAMPLING one it could have enumerated within budget.

That means the search half can be closed **by construction** rather than by accumulating
instruments: enumerate every latent configuration, take the EXACT argmin of the deployed
objective, and see whether it beats the pool.  If the exact optimum of the search does not beat a
zero-evaluation retrieval pool, no sampler -- quantum, classical, trained or not -- can, because
there is nothing left in the space to find.

THE DECODER IS STOCHASTIC, AND THAT DICTATES THE DESIGN.  `draw_from_basins` maps a bitstring to
a **von Mises draw**, not to a structure, so `2**n` bitstrings are `2**n` *distributions*.  The
deterministic, exhaustive object is therefore each bitstring decoded at its **von Mises MODE**
(`mu`), which is the latent's own mode set and the thing an argmin-readout VQE is choosing among.
A stochastic arm is reported beside it so the mode restriction is priced rather than assumed.

ARMS, per target.  Every arm native-free unless labelled ORACLE.

    latent_argmin    EXACT argmin of the deployed objective over all 2**n modes    the search optimum
    latent_oracle    the best of the same 2**n modes by TRUE RMSD                  ORACLE, the
                                                                                  latent's own ceiling
    latent_draw      one stochastic draw per bitstring, argmin by objective        prices the mode
                                                                                  restriction
    latent_mean      RMSD of the whole mode set                                    what the space
                                                                                  looks like on average
    pool_argmin      the shipped distogram score's argmin over the K=500 pool      the comparator
    pool_tail_avg    the shipped top-75 coordinate average                         THE INCUMBENT
                                                                                  (point cloud)

THE DECOMPOSITION THIS BUYS, and it is the point.  Two exhaustive numbers separate the two halves
of the problem cleanly, with no sampling error in either:

    latent_oracle  BAD   ->  the latent does not CONTAIN good structures.  Generation is the wall.
    latent_oracle  GOOD
      and argmin   BAD   ->  the space contains them and the objective cannot find them.
                            DISCRIMINATION is the wall, and it is proved by exhaustion.

PRE-REGISTRATION.

  Hypothesis:   the exact latent argmin does not beat the pool, and the latent's ORACLE ceiling
                is far better than its argmin -- i.e. discrimination binds, by exhaustion.
  Prediction:   latent_argmin >= pool_tail_avg (worse or equal); latent_oracle << latent_argmin.
  Primary:      latent_argmin - pool_tail_avg, paired over the enumerable targets, with SE and a
                per-comparison MDE = 2.8016 * SE reported beside it (L8: MDE is an operator).
  Falsifier:    if latent_argmin BEATS pool_tail_avg past that comparison's own MDE, the search
                half is NOT closed and a better sampler on this latent is worth building.
  Null:         `latent_mean` -- the whole space averaged.  If argmin does not beat the space's
                own mean, the objective has no ordering skill on the latent at all.
  Matched ctrl: `pool_argmin` at the same readout (single structure, argmin by score) -- so the
                latent and the pool are compared through the SAME operator, not one through an
                argmin and the other through an average.
  Budget:       exhaustive; 2**n <= 8192 evaluations per target, no more than the samplers got.
  Promotion:    none.  This is a bound.

  THE TRAP.  `latent_argmin` is a SINGLE STRUCTURE and `pool_tail_avg` is an AVERAGE of 75.  That
  is the same basis error the sprint has already made twice, so BOTH readouts are reported for
  BOTH sources and the primary is stated against the matched one.  Sprint 21 L3 measured the
  averaging operator at -0.28 to -0.31 on this instrument INDEPENDENT of the energy, so comparing
  an argmin to an average would hand the pool that margin for free.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)


def _save(obj, name=None):
    """ATOMIC.  Workstream D read a 45-row partial of this file mid-write and computed a
    statistic on it before noticing.  Every other lane's writer (d_enc._save, tailprice,
    a_matrix._write, qb2_run.ck_save) already does tmp + os.replace; this one did not."""
    path = os.path.join(RESULTS, name or ("latentfull%s.json" % TAG))
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)

from s12 import instrument as I             # noqa: E402
from s15 import seed as SD                  # noqa: E402
from s19 import qb_lib as QB                # noqa: E402

MIN_N = int(os.environ.get("LF_MIN_N", 0))
MAX_N = int(os.environ.get("LF_MAX_N", 13))   #: 2**13 = 8192, the budget the samplers received
TAG = os.environ.get("LF_TAG", "")            #: distinct output file for a distinct n-range
CHUNK = 2048


def run(targets=None):
    tg = targets if targets is not None else I.targets()
    tg = [t for t in tg if MIN_N <= int(t["n"]) <= MAX_N]
    rows, t0 = [], time.time()
    print(f"enumerable targets ({MIN_N} <= n <= {MAX_N}): {len(tg)}", flush=True)

    for c, t in enumerate(tg):
        pdb, n = t["pdb"], int(t["n"])
        tgt = QB.target(pdb)
        nat = tgt["nat"]
        mu = tgt["mu"]                                  # (n, 2, 2)
        N = 1 << n

        #: ---- EXHAUSTIVE: every bitstring, decoded at its von Mises MODE (deterministic)
        best_o, best_r, oracle_r, sum_r, cnt = np.inf, np.nan, np.inf, 0.0, 0
        rng = SD.stable_rng(pdb, "s21latentfull")
        best_draw_o, best_draw_r = np.inf, np.nan
        best_ob, best_rb = np.inf, np.nan
        #: `Obj.raw` IS the deployed objective (s19/qb_lib.py:236), called unbudgeted for
        #: post-hoc scoring exactly as its docstring provides.  Not a reimplementation.
        ob = QB.new_obj(tgt)
        for s0 in range(0, N, CHUNK):
            k = np.arange(s0, min(s0 + CHUNK, N), dtype=np.int64)
            bits = ((k[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)
            phi = mu[np.arange(n)[None, :], bits, 0]
            psi = mu[np.arange(n)[None, :], bits, 1]
            e, CA = ob.raw(phi, psi)
            #: the DEPLOYED selector, on the same structures.  Workstream D measured
            #: argmin_sq - argmin_bayes = +0.050 [-0.019,+0.156] (MDE 0.124, median rho 0.973):
            #: a null, but a KNOWN-DIRECTION null, so scoring the latent by the squared
            #: functional while the pool is scored by the Bayes risk would be an unstated
            #: operator difference INSIDE the primary comparison.  Both are carried.
            eb = np.asarray(I.shipped_score(tgt["dg"], I.pair_dists(CA, tgt["i"], tgt["j"])), float)
            rr = I.kabsch_rmsd_batch(CA, nat)           # ORACLE, scoring only
            a = int(np.argmin(e))
            if e[a] < best_o:
                best_o, best_r = float(e[a]), float(rr[a])
            ab = int(np.argmin(eb))
            if eb[ab] < best_ob:
                best_ob, best_rb = float(eb[ab]), float(rr[ab])
            oracle_r = min(oracle_r, float(rr.min()))
            sum_r += float(rr.sum()); cnt += len(k)
            #: stochastic arm -- one draw per bitstring, prices the mode restriction
            dphi, dpsi = QB.draw_from_basins(bits, mu, tgt["kap"], rng)
            ed, CAd = ob.raw(dphi, dpsi)
            rrd = I.kabsch_rmsd_batch(CAd, nat)
            ad = int(np.argmin(ed))
            if ed[ad] < best_draw_o:
                best_draw_o, best_draw_r = float(ed[ad]), float(rrd[ad])

        #: ---- the comparator, through BOTH readouts so the basis is matched
        u = I.load_univ(pdb)
        p = I.pool_idx(u)
        W = np.asarray(u["W"], float)[p]
        D = I.pair_dists(W, tgt["i"], tgt["j"])
        sc = np.asarray(I.shipped_score(tgt["dg"], D), float)
        order = np.argsort(sc, kind="stable")
        pool_argmin = float(I.ca_rmsd(W[order[0]], nat))
        avg, _b = I.coordinate_average(W[order[:75]])
        pool_tail_avg = float(I.ca_rmsd(np.asarray(avg, float), nat))
        #: BASIS-MATCHED comparator (Workstream D).  The latent arms are BUILT CHAINS -- torsions
        #: through build_ca_exact -- while `W` is each retrieved window's OWN coordinates.  That
        #: is a second unstated basis difference inside the primary, worth ~+0.011 A per member
        #: AGAINST the latent arm.  Here the pool's windows are put on the latent's own manifold:
        #: rebuilt from their torsions, then scored AND measured on the rebuild.
        Wr = np.asarray(I.build_ca(u["PHI"][p], u["PSI"][p]), float)
        scr = np.asarray(I.shipped_score(tgt["dg"], I.pair_dists(Wr, tgt["i"], tgt["j"])), float)
        pool_argmin_rb = float(I.ca_rmsd(Wr[int(np.argmin(scr))], nat))

        rows.append({"pdb": pdb, "n": n, "fold": int(t["fold"]), "N_latent": N,
                     "latent_argmin": best_r, "latent_argmin_bayes": best_rb,
                     "latent_oracle": oracle_r,
                     "latent_mean": sum_r / max(cnt, 1), "latent_draw": best_draw_r,
                     "pool_argmin": pool_argmin, "pool_argmin_rb": pool_argmin_rb,
                     "pool_tail_avg": pool_tail_avg})
        if (c + 1) % 5 == 0:
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            _save({"rows": rows, "complete": False, "max_n": MAX_N, "n_expected": len(tg)})

    ok = (len(rows) == len(tg)) and len(tg) > 0 and all(
        np.isfinite(r["latent_argmin"]) and np.isfinite(r["latent_argmin_bayes"])
        and np.isfinite(r["pool_argmin_rb"]) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "max_n": MAX_N, "n_expected": len(tg)})
    report(rows)
    return rows


def _stat(d, rng, B=4000):
    """mean, SE, per-comparison MDE (L8: MDE is an operator, not a constant), and a CI."""
    d = np.asarray(d, float); k = len(d)
    se = float(d.std(ddof=1) / np.sqrt(k))
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return (float(d.mean()), se, 2.8016 * se,
            float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5)),
            int((d < 0).sum()), int((d > 0).sum()))


def report(rows=None):
    if rows is None:
        rows = []
        for t in ([TAG] if TAG else ["", "_hi"]):
            f = os.path.join(RESULTS, "latentfull%s.json" % t)
            if os.path.exists(f):
                rows += json.load(open(f))["rows"]
    rng = SD.stable_rng("latentfull", "report")
    g = lambda k: np.array([r[k] for r in rows])           # noqa: E731
    print(f"\nn = {len(rows)} enumerable targets (n <= {MAX_N}).  EXHAUSTIVE over 2**n modes.")
    print(f"  latent sizes: min {int(g('N_latent').min())}  median "
          f"{int(np.median(g('N_latent')))}  max {int(g('N_latent').max())}\n")
    print(f"  {'arm':<16}{'basis':<14}{'RMSD':>8}{'median':>9}")
    for k, b in (("latent_mean", "mode set"), ("latent_argmin", "single"),
                 ("latent_argmin_bayes", "single"),
                 ("latent_draw", "single"), ("latent_oracle", "single ORACLE"),
                 ("pool_argmin", "single window"), ("pool_argmin_rb", "single REBUILT"),
                 ("pool_tail_avg", "avg of 75")):
        print(f"  {k:<16}{b:<14}{g(k).mean():>8.3f}{np.median(g(k)):>9.3f}")

    print("\n  PRIMARY, matched readout (single structure vs single structure):")
    for lab, a, b in (("latent_argmin_bayes - pool_argmin_rb  [MATCHED FUNC + BASIS]",
                       "latent_argmin_bayes", "pool_argmin_rb"),
                      ("latent_argmin_bayes - pool_argmin  [matched func, window basis]",
                       "latent_argmin_bayes", "pool_argmin"),
                      ("latent_argmin - pool_argmin  [squared vs Bayes]", "latent_argmin", "pool_argmin"),
                      ("latent_argmin - pool_tail_avg  [UNMATCHED, avg]", "latent_argmin", "pool_tail_avg")):
        m, se, mde, lo, hi, w, l = _stat(g(a) - g(b), rng)
        print(f"    {lab:<48}{m:+.3f}  SE {se:.3f}  MDE {mde:.3f}  "
              f"[{lo:+.3f},{hi:+.3f}]  {w}W/{l}L")

    print("\n  THE DECOMPOSITION -- does the latent CONTAIN good structures?")
    m, se, mde, lo, hi, w, l = _stat(g("latent_argmin") - g("latent_oracle"), rng)
    print(f"    argmin - ORACLE ceiling of the SAME 2**n modes {m:+.3f}  SE {se:.3f}  MDE {mde:.3f}")
    m2, se2, mde2, lo2, hi2, w2, l2 = _stat(g("latent_oracle") - g("pool_tail_avg"), rng)
    print(f"    latent ORACLE - incumbent                      {m2:+.3f}  SE {se2:.3f}  MDE {mde2:.3f}")
    m3, se3, mde3, *_ = _stat(g("latent_argmin") - g("latent_mean"), rng)
    print(f"    argmin - the space's own MEAN (the null)       {m3:+.3f}  SE {se3:.3f}  MDE {mde3:.3f}")
    m4, se4, mde4, *_ = _stat(g("latent_draw") - g("latent_argmin"), rng)
    print(f"    stochastic draw - mode argmin (mode price)     {m4:+.3f}  SE {se4:.3f}  MDE {mde4:.3f}")

    print("\nREAD.  If the latent's ORACLE ceiling is far below its argmin, the space CONTAINS good")
    print("structures the objective cannot find -- discrimination binds, proved BY EXHAUSTION rather")
    print("than by instrument count.  If the ceiling is also poor, generation is the wall instead.")
    print("The primary is the MATCHED readout; the unmatched row is printed only to show the")
    print("averaging margin it would otherwise hand the pool for free.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report()
    else:
        run()
