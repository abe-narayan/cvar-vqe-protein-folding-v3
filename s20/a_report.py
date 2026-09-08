"""SPRINT 20, AGENT A -- the report.  Reads `a_src.json` (point cloud) and `a_proj.json` (built).

Every table states its BASIS.  Nothing here compares a built chain to the 3.048 A point cloud.

Run:  python -m s20.a_report
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I              # noqa: E402
from s20 import a_src as S                   # noqa: E402

SRC = ["pool", "pep", "prot", "halfA", "halfB", "rand500", "tors", "unsel", "helix"]


def _load(name):
    """Read a result file that a live checkpointing run may be rewriting under us.

    File contention corrupted a Sprint-19 lane's read; a partial parse must fail loudly and
    retry, never be silently treated as a short run.
    """
    p = os.path.join(S.RESULTS, name)
    if not os.path.exists(p):
        return None, None
    for _ in range(12):
        try:
            d = json.load(open(p))
            return d["rows"], bool(d.get("complete"))
        except (json.JSONDecodeError, KeyError):
            time.sleep(1.0)
    raise RuntimeError(f"{name} unreadable after 12 attempts (concurrent writer?)")


def rho(rows, a, b, key="rho"):
    for x, y in ((a, b), (b, a)):
        k = f"{key}|{x}|{y}"
        if k in rows[0]:
            return np.array([r[k] for r in rows], float)
    return None


def main():
    rows, comp = _load("a_src.json")
    n = len(rows)
    folds = np.array([r["fold"] for r in rows], int)
    pdbs = [r["pdb"] for r in rows]
    g = lambda k: np.array([r[k] for r in rows], float)      # noqa: E731
    print(f"=== S20 AGENT A -- n = {n}  complete = {comp} ===")
    print(f"top-75 reconstruction set-identical to shipped_record['sub']: "
          f"{json.load(open(os.path.join(S.RESULTS,'a_src.json')))['gate_ok']}/{n}")
    print(f"peptide fraction of the shipped 500-pool: {g('pep_frac_univ').mean():.3f} of the "
          f"universe, {g('pep_frac_pool').mean():.3f} of the pool\n")

    # ---------------------------------------------------------------- sources
    print("--- SOURCES.  BASIS: POINT CLOUD (coordinate average).  reference pool = 3.048 A ---")
    hd = (f"{'source':<9}{'realis':>8}{'gen ceil':>9}{'sel ceil':>9}{'member':>8}"
          f"{'cov<2A':>8}{'divGEO':>8}{'divERR':>8}{'bond':>7}{'bondmin':>8}"
          f"{'rho|pool':>10}{'memrho':>8}")
    print(hd)
    tab = {}
    for s in SRC:
        r_ = rho(rows, s, "pool")
        mr = f"{g('memrho|pool|'+s).mean():>8.3f}" if f"memrho|pool|{s}" in rows[0] else f"{'-':>8}"
        if s == "helix":
            print(f"{s:<9}{g('rmsd|'+s).mean():>8.3f}{'-':>9}{'-':>9}{'-':>8}{'-':>8}"
                  f"{'-':>8}{'-':>8}{g('bond|'+s).mean():>7.3f}{g('bondmin|'+s).mean():>8.3f}"
                  f"{np.nanmean(r_):>10.3f}{'-':>8}")
            tab[s] = {"rmsd_cloud": float(g("rmsd|" + s).mean()),
                      "rho_pool": float(np.nanmean(r_))}
            continue
        print(f"{s:<9}{g('rmsd|'+s).mean():>8.3f}{g('gen_ceiling|'+s).mean():>9.3f}"
              f"{g('sel_ceiling|'+s).mean():>9.3f}{g('member_rmsd|'+s).mean():>8.3f}"
              f"{g('cov2|'+s).mean():>8.3f}{g('div_geo|'+s).mean():>8.3f}"
              f"{g('div_err|'+s).mean():>8.3f}{g('bond|'+s).mean():>7.3f}"
              f"{g('bondmin|'+s).mean():>8.3f}"
              f"{(np.nanmean(r_) if r_ is not None else np.nan):>10.3f}{mr}")
        tab[s] = {"rmsd_cloud": float(g("rmsd|" + s).mean()),
                  "gen_ceiling": float(g("gen_ceiling|" + s).mean()),
                  "sel_ceiling": float(g("sel_ceiling|" + s).mean()),
                  "member_rmsd": float(g("member_rmsd|" + s).mean()),
                  "cov_under2": float(g("cov2|" + s).mean()),
                  "div_geo": float(g("div_geo|" + s).mean()),
                  "div_err": float(g("div_err|" + s).mean()),
                  "bond": float(g("bond|" + s).mean()),
                  "bond_min": float(g("bondmin|" + s).mean()),
                  "bond_member": float(g("bondmem|" + s).mean()),
                  "rho_pool": float(np.nanmean(r_)) if r_ is not None else None}
    print(f"{'NATIVE':<9}{0.0:>8.3f}{'':>9}{'':>9}{'':>8}{'':>8}{'':>8}{'':>8}"
          f"{g('bond|native').mean():>7.3f}")

    # ---------------------------------------------------------------- the falsifier
    print("\n--- F-A1: IS THE CORPUS THE CARRIER? ---")
    pp = rho(rows, "pep", "prot")
    hh = rho(rows, "halfA", "halfB")
    pep_pool = rho(rows, "pep", "pool")
    ppL, hhL = rho(rows, "pep", "prot", "rhoL"), rho(rows, "halfA", "halfB", "rhoL")
    print(f"  rho(e_pep, e_pool)              = {np.nanmean(pep_pool):.4f}   "
          f"median {np.nanmedian(pep_pool):.4f}   [F-A1a bar: >= 0.85]")
    print(f"  rho(e_pep, e_prot)  CORPUS split= {np.nanmean(pp):.4f}   median {np.nanmedian(pp):.4f}")
    print(f"  rho(e_halfA,e_halfB) RANDOM split= {np.nanmean(hh):.4f}   median "
          f"{np.nanmedian(hh):.4f}   <- the matched ceiling")
    print(f"  share of ceiling = {np.nanmean(pp)/np.nanmean(hh):.4f}   [F-A1b bar: >= 0.85]")
    st = I.paired(pp, hh, folds=folds, names=pdbs)
    print(f"  corpus - random  = {st['mean_diff']:+.4f} [{st['ci95'][0]:+.4f},{st['ci95'][1]:+.4f}]"
          f"  {st['n_better']}W/{st['n_worse']}L (W = corpus decorrelates more)")
    print(f"  long-range only (sep>=5): corpus {np.nanmean(ppL):.4f}  random {np.nanmean(hhL):.4f}")
    fa1a = float(np.nanmean(pep_pool)) >= 0.85
    fa1b = float(np.nanmean(pp) / np.nanmean(hh)) >= 0.85
    print(f"  >>> F-A1a {'FIRES' if fa1a else 'does not fire'}   "
          f"F-A1b {'FIRES' if fa1b else 'does not fire'}")

    print("\n--- cross-source structural agreement, BASIS: POINT CLOUD ---")
    print(f"  {'pair':<22}{'RMSD(X_A,X_B)':>15}{'each to native':>16}")
    for a, b in [("pool", "pep"), ("pool", "prot"), ("pool", "halfA"), ("halfA", "halfB"),
                 ("pep", "prot"), ("pool", "rand500"), ("pool", "tors"), ("pool", "unsel"),
                 ("pool", "helix")]:
        k = f"xrmsd|{a}|{b}" if f"xrmsd|{a}|{b}" in rows[0] else f"xrmsd|{b}|{a}"
        print(f"  {a+'~'+b:<22}{g(k).mean():>15.3f}"
              f"{g('rmsd|'+a).mean():>9.3f}/{g('rmsd|'+b).mean():.3f}")

    # ---------------------------------------------------------------- fusion
    print("\n--- FUSION, BASIS: POINT CLOUD.  vs the `pool` arm (3.048 A reference) ---")
    base = g("rmsd|pool")
    for k in sorted([k for k in rows[0] if k.startswith(("fuse|", "merge|"))]):
        v = g(k)
        st = I.paired(v, base, folds=folds, names=pdbs)
        bk = ("fusebond|" + k.split("|", 1)[1]) if k.startswith("fuse|") \
            else ("mergebond|" + k.split("|", 1)[1])
        print(f"  {k:<24}{v.mean():>7.3f}  {st['mean_diff']:+.3f} "
              f"[{st['ci95'][0]:+.3f},{st['ci95'][1]:+.3f}]  "
              f"{st['n_better']}W/{st['n_worse']}L  bond {g(bk).mean():.3f}")

    # ---------------------------------------------------------------- built chains
    prows, pcomp = _load("a_proj.json")
    if prows:
        pf = np.array([r["fold"] for r in prows], int)
        pp_ = [r["pdb"] for r in prows]
        h = lambda k: np.array([r[k] for r in prows], float)   # noqa: E731
        print(f"\n--- BUILT CHAINS (lam=0 projection).  BASIS: BUILT.  incumbent 3.204 A.  "
              f"n = {len(prows)} complete={pcomp} ---")
        print(f"  {'arm':<24}{'cloud':>8}{'BUILT':>8}{'tax':>7}{'bond':>7}"
              f"{'   vs pool(built)'}")
        keys = [k[6:] for k in prows[0] if k.startswith("rmsdP|")]
        bb = h("rmsdP|pool")
        for k in keys:
            v = h("rmsdP|" + k)
            st = I.paired(v, bb, folds=pf, names=pp_)
            print(f"  {k:<24}{h('rmsdC|'+k).mean():>8.3f}{v.mean():>8.3f}"
                  f"{(v-h('rmsdC|'+k)).mean():>7.3f}{h('bondP|'+k).mean():>7.3f}"
                  f"   {st['mean_diff']:+.3f} [{st['ci95'][0]:+.3f},{st['ci95'][1]:+.3f}]"
                  f"  {st['n_better']}W/{st['n_worse']}L")
        tab["built"] = {k: {"cloud": float(h("rmsdC|" + k).mean()),
                            "built": float(h("rmsdP|" + k).mean()),
                            "bond": float(h("bondP|" + k).mean())} for k in keys}
        tab["built_n"] = len(prows)
        tab["built_complete"] = pcomp

        print(f"\n--- THE FALSIFIER ON THE BUILT BASIS (every chain at bond 3.804 A, so no "
              f"alignment here can be a contraction artefact).  n = {len(prows)} ---")
        print(f"  {'arm':<24}{'rhoB vs pool':>14}{'partialled on helix':>22}"
              f"{'RMSD(X,X_pool)':>16}")
        for k in ["pep", "prot", "halfA", "halfB", "rand500", "tors", "unsel", "helix"]:
            if "rhoB|" + k not in prows[0]:
                continue
            rp = f"{h('rhoBP|'+k).mean():>22.3f}" if "rhoBP|" + k in prows[0] else f"{'-':>22}"
            xr = f"{h('xrmsdB|'+k).mean():>16.3f}" if "xrmsdB|" + k in prows[0] else f"{'-':>16}"
            print(f"  {k:<24}{h('rhoB|'+k).mean():>14.3f}{rp}{xr}")
        for k, lab in [("pep|prot", "CORPUS partition"),
                       ("halfA|halfB", "RANDOM partition (matched ceiling)")]:
            if "rhoB|" + k in prows[0]:
                print(f"  {k:<24}{h('rhoB|'+k).mean():>14.3f}"
                      f"{h('rhoBP|'+k).mean():>22.3f}{'':>16}   {lab}")
        if "rhoB|pep|prot" in prows[0]:
            sh = h("rhoB|pep|prot").mean() / h("rhoB|halfA|halfB").mean()
            shp = h("rhoBP|pep|prot").mean() / h("rhoBP|halfA|halfB").mean()
            print(f"  >>> BUILT-basis share of ceiling: raw {sh:.4f}  partialled {shp:.4f}"
                  f"   [F-A1b bar 0.85]")
            tab["falsifier_built"] = {"rho_pep_pool": float(h("rhoB|pep").mean()),
                                      "rho_pep_prot": float(h("rhoB|pep|prot").mean()),
                                      "rho_halfA_halfB": float(h("rhoB|halfA|halfB").mean()),
                                      "share_raw": float(sh), "share_partial": float(shp),
                                      "n": len(prows), "complete": pcomp}

    tab["falsifier"] = {"rho_pep_pool": float(np.nanmean(pep_pool)),
                        "rho_pep_prot": float(np.nanmean(pp)),
                        "rho_halfA_halfB": float(np.nanmean(hh)),
                        "share_of_ceiling": float(np.nanmean(pp) / np.nanmean(hh)),
                        "F_A1a_fires": bool(fa1a), "F_A1b_fires": bool(fa1b)}
    json.dump(tab, open(os.path.join(S.RESULTS, "a_report.json"), "w"), indent=1, default=float)
    return tab


if __name__ == "__main__":
    main()
