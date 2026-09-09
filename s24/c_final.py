"""s24/c_final.py -- LANE C ADDENDUM: recipe reconciliation, borrowed audits, provenance stamps.

Three things the coordinator asked for after `c_ladder.py` had already run, done in one pass so the
promoted artefacts are replicable from their own hashes.

1. RECIPE RECONCILIATION FOR LANE D.  D computes the same distogram-alignment quantity on its own
   30-target panel using `s24/referent.py::_stats` verbatim.  Four of my five choices match D's.
   ONE DIFFERS AND IT IS THE FIRST ONE:

       (a) Dhat      D uses the MAP, grid[argmin(risk, axis=1)].  `c_ladder.py` used
                     dg["expected"], the POSTERIOR MEAN.  THEY ARE DIFFERENT ESTIMATORS.
       (b) min_sep   both 2 (the instrument default, via I.pair_index)                    MATCH
       (c) Dc        both the uniform coordinate average in its own medoid frame, point
                     cloud -- not the medoid member, not projected                        MATCH
       (d) averaging both per target then averaged, never pooled over all pairs first     MATCH
       (e) centring  neither eC nor eP centred                                            MATCH

   So this file recomputes the alignment under BOTH estimators, on all 126, so D can diff against
   my per-target numbers on any of its 30 without either of us guessing.  Reporting only my own
   choice and letting D assume it matched is exactly the failure mode that makes two lanes'
   numbers look like a replication when they are two different measurements.

2. AUDITS FROM LANE B'S HARNESS RATHER THAN MINE.  `s24/residlib.py` is complete and exercised, so
   `collapse_audit` and `validity_audit` are taken from it verbatim.  Per Lane B's finding, for
   structures built from torsions omega deviation, cis fraction, bond length and bond angle
   deviation are ZERO BY CONSTRUCTION and Ca-Ca is 3.80 A by construction; those axes are named
   and omitted rather than printed as a table of structural zeros.  **The only validity axes
   carrying information in torsion space are RAMACHANDRAN and CLASHES**, and those are reported.

   The emitted sets are REPLAYED from the identical seed stream `c_ladder.py` used, so this audits
   the sets whose RMSD is quoted, not a fresh draw from the same sampler.

3. PROVENANCE.  Every artefact this lane promotes is re-saved through
   `stats_lib.save_atomic(..., module_file=...)`, which stamps module name, source sha256, git
   commit, dirty flag and launch time, and sets `complete` on the FULL KEY SET rather than a row
   count.  `c_ladder.json` is stamped with the hash of `c_ladder.py`, which has not been edited
   since that run.

ORACLE labelling unchanged: every RMSD and every alignment cosine reads the native and is post-hoc
scoring of a native-free decision.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402
from s24 import stats_lib as ST          # noqa: E402
from s24 import residlib as RL           # noqa: E402
from s24 import c_probe as CP            # noqa: E402
from s24 import c_ladder as CL           # noqa: E402
from core import project as pj           # noqa: E402

TOPM = CL.TOPM
NSAMP = CL.NSAMP
ARMS = CL.ARMS
NAUD = 250                    # subsample for the O(B^2) collapse audit; stated, not hidden
OUT = os.path.join(RES, "c_final.json")
NEED = ("pdb", "map", "exp", "audit")


def dhat_map(dg):
    """Lane D's recipe (a): the MAP of the per-pair risk grid."""
    return np.asarray(dg["grid"], float)[np.argmin(np.asarray(dg["risk"], float), axis=1)]


def native_ss(nat):
    """CA-ONLY secondary-structure proxy for a NATIVE trace.  Stated as a proxy, not DSSP.

    Lane A's L8 found the permitted corpus is H 0.375 / E 0.022 / C 0.604 -- beta sheet is
    essentially ABSENT, and mechanically so: DSSP assigns E only with a paired strand, and both
    banks hold isolated short chains whose partner strand was left outside the window by
    construction.  So every arm in this ladder is helix/coil-biased and must be expected to fail
    quietly on sheet targets.  Stratifying by native SS is what makes that visible instead of
    averaged away.

    Only `nat_ca` is available for the dev natives (no backbone N/C), so DSSP is not computable
    here.  The standard CA-only discriminator is used instead: the i->i+3 and i->i+4 spans, which
    separate a 3.6-residue helical turn (tight) from an extended strand (long).
    """
    nat = np.asarray(nat, float); n = len(nat)
    if n < 5:
        return "coil", 0.0, 0.0
    d3 = np.linalg.norm(nat[3:] - nat[:-3], axis=-1)          # i -> i+3
    d4 = np.linalg.norm(nat[4:] - nat[:-4], axis=-1)          # i -> i+4
    m = len(d4)
    helical = (d3[:m] > 4.5) & (d3[:m] < 6.5) & (d4 < 7.2)
    h = float(helical.mean())
    e = float((d3 > 9.0).mean())
    cls = "helix" if h >= 0.4 else ("extended" if e >= 0.4 else "coil")
    return cls, h, e


def run():
    tg = I.targets()
    rows = []
    for ci, t in enumerate(tg):
        pdb = t["pdb"]; n = t["n"]
        u = I.load_univ(pdb); nat = u["nat_ca"]
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(n)      # min_sep=2 (b)
        dnat = I.pair_dists(nat[None], i, j)[0]
        g_map = dhat_map(dg) - dnat                       # (a) D's estimator
        g_exp = np.asarray(dg["expected"], float) - dnat  # (a) mine

        idx = I.pool_idx(u); Wp = u["W"][idx]
        scp = np.asarray(I.shipped_score(dg, I.pair_dists(Wp, i, j)), float)
        top = np.argsort(scp, kind="stable")[:TOPM]
        A = Wp[top]
        cA = CP._avg(A)                                   # (c) uniform average, own medoid frame
        eA = I.pair_dists(cA[None], i, j)[0] - dnat       # (e) uncentred

        ss, hfrac, efrac = native_ss(nat)
        r = {"pdb": pdb, "n": n, "fold": int(u["fold"]),
             "ss": ss, "h_frac": hfrac, "e_frac": efrac,
             "map": {"INCUMBENT": CP._cos(eA, g_map)},
             "exp": {"INCUMBENT": CP._cos(eA, g_exp)},
             "audit": {}}

        # --- replay c_ladder.py's EXACT seed stream so this audits the quoted sets
        rng0 = SD.stable_rng("c_ladder", pdb)
        for tag in ARMS:
            def fn(g, _tag=tag):
                return {"T0_helix": lambda: CP.s_helix(n, NSAMP, g),
                        "T1_blind": lambda: CP.s_blind(u, n, NSAMP, g),
                        "T2_restype": lambda: CP.s_restype(u, n, NSAMP, g),
                        "T3_pool": lambda: CP.s_pool(u, n, NSAMP, g, idx[top])}[_tag]()

            ph, ps = fn(rng0)
            CA = np.asarray(pj.build_ca_exact(ph, ps), float)
            sc = np.asarray(I.shipped_score(dg, I.pair_dists(CA, i, j)), float)
            o = np.argsort(sc, kind="stable")
            cS = CP._avg(CA[o[:TOPM]])
            ru = rng0.choice(len(CA), TOPM, replace=False)
            cU = CP._avg(CA[ru])
            _ph2, _ps2 = fn(rng0)                          # consume draw 2 exactly as c_ladder did
            _CA2 = np.asarray(pj.build_ca_exact(_ph2, _ps2), float)
            _ru2 = rng0.choice(len(_CA2), TOPM, replace=False)

            eS = I.pair_dists(cS[None], i, j)[0] - dnat
            eU = I.pair_dists(cU[None], i, j)[0] - dnat
            r["map"][tag] = {"scored": CP._cos(eS, g_map), "unif": CP._cos(eU, g_map)}
            r["exp"][tag] = {"scored": CP._cos(eS, g_exp), "unif": CP._cos(eU, g_exp)}

            # --- Lane B's harness, on a stated subsample (collapse_audit is O(B^2))
            s = SD.stable_rng("c_final", pdb, tag).permutation(len(CA))[:NAUD]
            r["audit"][tag] = {"collapse": RL.collapse_audit(ph[s], ps[s], CA[s]),
                               "validity": RL.validity_audit(ph, ps)}
        rows.append(r)
        if (ci + 1) % 10 == 0:
            print("  %d/%d" % (ci + 1, len(tg)), flush=True)
            ST.save_atomic(OUT, {"rows": rows, "n_audit_subsample": NAUD},
                           complete_keys=NEED, rows=rows, n_expected=len(tg),
                           module_file=__file__)
    ST.save_atomic(OUT, {"rows": rows, "n_audit_subsample": NAUD, "nsamp": NSAMP,
                         "recipe": {"dhat": "BOTH map and expected reported", "min_sep": 2,
                                    "Dc": "uniform coordinate average, own medoid frame, point cloud",
                                    "averaging": "per target then averaged", "centring": "none"}},
                   complete_keys=NEED, rows=rows, n_expected=len(tg), module_file=__file__)
    stamp_existing()
    report(rows)
    return rows


def stamp_existing():
    """Re-save the lane's already-complete artefacts through save_atomic, stamping the source that
    actually produced them.  Contents are NOT recomputed -- only provenance and the completion flag
    on the full key set are added."""
    for name, mod, need in (("c_ladder.json", os.path.join(HERE, "c_ladder.py"),
                             ("pdb", "incumbent", "rebuild75") + ARMS),
                            ("c_bias.json", os.path.join(HERE, "c_bias.py"),
                             ("pdb", "dgram_signed", "members_signed", "average_signed"))):
        p = os.path.join(RES, name)
        if not os.path.exists(p):
            continue
        o = json.load(open(p)); rows = o.get("rows", [])
        ST.save_atomic(p, o, complete_keys=need, rows=rows, n_expected=126, module_file=mod)
        print("  stamped %s  complete=%s" % (name, json.load(open(p))["complete"]))


def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    print("\n" + "=" * 92)
    print("LANE C ADDENDUM, n=%d." % len(rows))
    print("=" * 92)

    print("\n[A] RECIPE RECONCILIATION.  cos(emitted distance error, distogram's own error),")
    print("    under Lane D's MAP estimator and under my posterior-mean estimator.")
    print("    Choices (b) min_sep=2, (c) uniform average in own medoid frame / point cloud,")
    print("    (d) per-target then averaged, (e) uncentred -- ALL MATCH D.  Only (a) differed.")
    print("\n  %-14s %12s %12s %10s | %12s %12s"
          % ("arm", "MAP scored", "EXP scored", "delta", "MAP unif", "EXP unif"))
    im = np.array([r["map"]["INCUMBENT"] for r in rows], float)
    ie = np.array([r["exp"]["INCUMBENT"] for r in rows], float)
    print("  %-14s %12.4f %12.4f %10.4f" % ("INCUMBENT", im.mean(), ie.mean(),
                                            im.mean() - ie.mean()))
    for k in ARMS:
        ms = np.array([r["map"][k]["scored"] for r in rows], float)
        es = np.array([r["exp"][k]["scored"] for r in rows], float)
        mu = np.array([r["map"][k]["unif"] for r in rows], float)
        eu = np.array([r["exp"][k]["unif"] for r in rows], float)
        print("  %-14s %12.4f %12.4f %10.4f | %12.4f %12.4f"
              % (k, ms.mean(), es.mean(), ms.mean() - es.mean(), mu.mean(), eu.mean()))
    dm = np.array([[r["map"][k]["scored"] - r["map"][k]["unif"] for k in ARMS] for r in rows])
    de = np.array([[r["exp"][k]["scored"] - r["exp"][k]["unif"] for k in ARMS] for r in rows])
    print("\n  THE LANE'S HEADLINE QUANTITY -- what SELECTION adds, averaged over the four arms:")
    print("    under the MAP            %+.4f" % dm.mean())
    print("    under the posterior mean %+.4f" % de.mean())
    print("    per arm (MAP): " + "  ".join("%s %+.3f" % (k[:2], dm[:, q].mean())
                                            for q, k in enumerate(ARMS)))

    print("\n[B] MODE-COLLAPSE AUDIT, Lane B's `residlib.collapse_audit`, %d-sample subsample."
          % NAUD)
    print("    Lane B reference values from 75-member sets: unique 66.8/75, dup 0.110, ESS 63.7,")
    print("    max-mode 0.044, pairRMSD 2.528, torsion entropy 0.261.")
    print("  %-14s %10s %9s %9s %10s %10s %9s"
          % ("arm", "uniq/250", "dup_frac", "ESS", "maxmode", "pairRMSD", "tors_H"))
    for k in ARMS:
        f = lambda q: np.array([r["audit"][k]["collapse"][q] for r in rows], float)   # noqa: E731
        print("  %-14s %10.1f %9.4f %9.1f %10.4f %10.4f %9.4f"
              % (k, f("n_unique").mean(), f("dup_frac").mean(), f("ess").mean(),
                 f("mode_occ_max").mean(), f("pair_rmsd_mean").mean(), f("tors_entropy").mean()))

    print("\n[C] GEOMETRIC VALIDITY.  Only RAMACHANDRAN and CLASHES carry information here:")
    print("    omega deviation, cis fraction, bond length and bond angle deviation are ZERO BY")
    print("    CONSTRUCTION for torsion-built chains, and Ca-Ca is 3.80 A by construction.")
    print("    Named and omitted rather than printed as a table of structural zeros.")
    print("  %-14s %14s %14s %14s" % ("arm", "rama_favoured", "clash_frac", "clash/struct"))
    for k in ARMS:
        f = lambda q: np.array([r["audit"][k]["validity"][q] for r in rows], float)   # noqa: E731
        print("  %-14s %14.4f %14.4f %14.4f"
              % (k, f("rama_favoured").mean(), f("clash_frac").mean(),
                 f("clash_per_struct").mean()))
    print("  Lane B's warning: a damped or mean-like sampler degrades EXACTLY this column while")
    print("  its RMSD column looks safest.  T0_helix (a jittered constant) is this lane's damped")
    print("  arm and is the one to read here.")

    ss_report(rows)


def ss_report(rows):
    """Stratify the confirmatory ladder by native SS class (Lane A L8: the corpus is H 0.375 /
    E 0.022 / C 0.604, so sheet targets are the predicted quiet failure)."""
    lad = os.path.join(RES, "c_ladder.json")
    if not os.path.exists(lad):
        return
    L = {r["pdb"]: r for r in json.load(open(lad))["rows"]}
    ss = {r["pdb"]: r["ss"] for r in rows}
    print("\n[D] STRATIFIED BY NATIVE SS CLASS (CA-only proxy, stated as a proxy, not DSSP).")
    print("    Lane A L8: the permitted corpus is H 0.375 / E 0.022 / C 0.604 -- beta sheet is")
    print("    essentially absent, so a sheet target is where this ladder should fail quietly.")
    classes = ["helix", "extended", "coil"]
    print("  %-10s %5s %10s %10s | %s" % ("class", "n", "incumbent", "poolbest",
                                          "  ".join("%-11s" % k for k in ARMS)))
    for c in classes:
        p = [q for q in L if ss.get(q) == c]
        if not p:
            continue
        inc = np.array([L[q]["incumbent"] for q in p], float)
        pb = np.array([L[q]["pool_best_ORACLE"] for q in p], float)
        cells = []
        for k in ARMS:
            v = np.array([L[q][k]["rmsd"] for q in p], float)
            cells.append("%6.3f%+5.2f" % (v.mean(), v.mean() - inc.mean()))
        print("  %-10s %5d %10.4f %10.4f | %s" % (c, len(p), inc.mean(), pb.mean(),
                                                  "  ".join(cells)))
    print("  Cells are the arm's standalone mean and its gap to the incumbent WITHIN that class.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
