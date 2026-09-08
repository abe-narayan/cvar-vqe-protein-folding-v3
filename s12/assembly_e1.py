"""E1 -- ORACLE floors of multi-piece fragment assembly (torsion concatenation, ideal geometry).

Everything in this file is ORACLE/DIAGNOSTIC: the native CA trace picks the pieces.
Per target we record:
  native_rebuild  : native torsions -> build_ca -> RMSD to native   (ideal-geometry floor)
  univ_best_real  : whole-window universe best on REAL coords (the 1.313 number)
  univ_best_ideal : same windows rebuilt from their torsions (m=1 control on the same geometry)
  local[L]        : mean over positions of the best local piece RMSD (ideal & real CA), L=3..10
  m=2 : for every cut c: per-piece argmin concat (E_concat), local errors, combination oracle
        (beam k=50 on global RMSD), junction-optimised placement floor (from argmin pieces),
        overlap-1 conventions A/B/avg (argmin and k=30 combination)
  m=3 : all compositions: argmin concat + combination k=15; junction floor on the equal split
  m=4 : equal split: argmin concat + combination k=8 + junction floor
"""
from __future__ import annotations
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import assembly_bank as AB
from s12 import assembly_common as AC

TOPK = 100


def local_table(fold, n, nat):
    """ORACLE: for every (s, L) the top-TOPK bank pieces by ideal-CA local RMSD, plus real-CA best."""
    tab = {}
    for L in range(AC.LMIN, AC.LMAX + 1):
        for s in range(0, n - L + 1):
            r = AC.local_rmsd_all(fold, s, L, nat[s:s + L])
            top = np.argsort(r, kind="stable")[:TOPK]
            rreal = AC.local_rmsd_all(fold, s, L, nat[s:s + L], real=True)
            tab[(s, L)] = (top, r[top], float(rreal.min()))
    return tab


def combos_global(fold, n, nat, comp, tab, k, overlap=0, conv="cut"):
    """Combination oracle: all k^m combinations of the top-k local pieces, global RMSD min.
    Returns (best_rmsd, best_choice(list of bank idx), argmin_concat_rmsd, mean_local_err)."""
    iv = AC.intervals(comp)
    tops, tors = [], []
    for q, (s, L) in enumerate(iv):
        if overlap and q > 0:
            s2, L2 = s - 1, L + 1
        else:
            s2, L2 = s, L
        top, rl, _ = tab[(s2, L2)]
        top = top[:k]
        ph, ps, _ = AB.piece_torsions(fold, L2, top)
        tops.append((s, L, top, rl[:k])); tors.append((ph, ps))
    m = len(iv)
    grids = np.array(np.meshgrid(*[np.arange(k)] * m, indexing="ij")).reshape(m, -1).T
    PHI = np.zeros((len(grids), n)); PSI = np.zeros((len(grids), n))
    for row, g in enumerate(grids):
        placed = [(tops[q][0], tops[q][1], tors[q][0][g[q]], tors[q][1][g[q]]) for q in range(m)]
        ph, ps = AC.concat_torsions(n, placed, overlap=overlap, conv=conv)
        PHI[row] = ph; PSI[row] = ps
    r = I.kabsch_rmsd_batch(AC.build_many(PHI, PSI), nat)
    b = int(np.argmin(r))
    argmin_r = float(r[0])           # row 0 = all per-piece argmins
    loc = float(np.mean([tops[q][3][0] for q in range(m)]))
    return float(r[b]), [int(tops[q][2][grids[b][q]]) for q in range(m)], argmin_r, loc, (PHI[0], PSI[0])


def run_target(t):
    pdb, n, fold = t["pdb"], t["n"], t["fold"]
    out_path = os.path.join(I.CACHE, f"asm_e1_{pdb}.json")
    if os.path.exists(out_path):
        return json.load(open(out_path))
    t0 = time.time()
    import peptide_db as db
    u = I.load_univ(pdb); nat = u["nat_ca"]
    q = db.by_pdb(pdb)
    res = {"pdb": pdb, "n": n, "fold": fold}
    res["native_rebuild"] = I.ca_rmsd(I.build_ca(q.phi[None], q.psi[None])[0], nat)
    res["univ_best_real"] = float(u["rr"].min())
    ru = I.kabsch_rmsd_batch(AC.build_many(u["PHI"], u["PSI"]), nat)
    res["univ_best_ideal"] = float(ru.min())
    res["pool_best_real"] = float(u["rr"][I.pool_idx(u)].min())
    tab = local_table(fold, n, nat)
    res["local"] = {L: {"ideal": float(np.mean([tab[(s, L)][1][0] for s in range(n - L + 1)])),
                        "real": float(np.mean([tab[(s, L)][2] for s in range(n - L + 1)])),
                        "ideal_top10": float(np.mean([tab[(s, L)][1][9] for s in range(n - L + 1)]))}
                    for L in range(AC.LMIN, AC.LMAX + 1)}
    # ---- m = 2
    m2 = {}
    for c in AC.cuts_2(n):
        comp = (c, n - c)
        best, choice, argmin_r, loc, (ph0, ps0) = combos_global(fold, n, nat, comp, tab, 50)
        d = {"argmin_concat": argmin_r, "combo50": best, "local_err": loc}
        for conv in ("A", "B", "avg"):
            if n - c + 1 <= AC.LMAX:
                b2, _, a2, l2, _ = combos_global(fold, n, nat, comp, tab, 30, overlap=1, conv=conv)
            else:
                b2 = a2 = float("nan")
            d[f"ov1_{conv}_argmin"] = a2; d[f"ov1_{conv}_combo30"] = b2
        m2[c] = d
    res["m2"] = m2
    cs = AC.cuts_2(n); mid = cs[len(cs) // 2] if len(cs) % 2 else cs[len(cs) // 2 - 1]
    bestc = min(cs, key=lambda c: m2[c]["combo50"])
    res["m2_mid_cut"] = mid; res["m2_best_cut"] = bestc
    for tag, c in (("mid", mid), ("best", bestc)):
        comp = (c, n - c)
        _, _, _, _, (ph0, ps0) = combos_global(fold, n, nat, comp, tab, 1)
        jr, _, _ = AC.optimise_junctions(ph0, ps0, comp, nat)
        res[f"m2_{tag}_junction_floor"] = jr
    # ---- m = 3
    m3 = {}
    for comp in AC.cuts_m(n, 3, 3, 8):
        best, choice, argmin_r, loc, _ = combos_global(fold, n, nat, comp, tab, 15)
        m3["-".join(map(str, comp))] = {"argmin_concat": argmin_r, "combo15": best, "local_err": loc}
    res["m3"] = m3
    eq3 = AC.equal_split(n, 3)
    _, _, _, _, (ph0, ps0) = combos_global(fold, n, nat, eq3, tab, 1)
    res["m3_equal"] = "-".join(map(str, eq3))
    res["m3_equal_junction_floor"], _, _ = AC.optimise_junctions(ph0, ps0, eq3, nat)
    # ---- m = 4 (equal split only)
    eq4 = AC.equal_split(n, 4)
    if min(eq4) >= AC.LMIN:
        best, choice, argmin_r, loc, (ph0, ps0) = combos_global(fold, n, nat, eq4, tab, 8)
        res["m4_equal"] = {"comp": "-".join(map(str, eq4)), "argmin_concat": argmin_r, "combo8": best, "local_err": loc}
        res["m4_equal"]["junction_floor"], _, _ = AC.optimise_junctions(ph0, ps0, eq4, nat)
    res["seconds"] = time.time() - t0
    json.dump(res, open(out_path, "w"), default=float)
    return json.load(open(out_path))


def _init():
    os.environ["OMP_NUM_THREADS"] = "2"
    try:
        import torch; torch.set_num_threads(2)
    except Exception:
        pass


def aggregate(rows):
    pdbs = [r["pdb"] for r in rows]
    G = lambda v: AC.group_means(v, pdbs)
    agg = {"n_targets": len(rows)}
    for key in ("native_rebuild", "univ_best_real", "univ_best_ideal", "pool_best_real",
                "m2_mid_junction_floor", "m2_best_junction_floor", "m3_equal_junction_floor"):
        agg[key] = G([r[key] for r in rows])
    agg["local"] = {L: {k: G([r["local"][str(L)][k] for r in rows]) for k in ("ideal", "real", "ideal_top10")}
                    for L in range(AC.LMIN, AC.LMAX + 1)}
    m2keys = ["argmin_concat", "combo50", "local_err"] + [f"ov1_{c}_{w}" for c in ("A", "B", "avg") for w in ("argmin", "combo30")]
    def GN(v):   # nan-aware group means (overlap arms are undefined when the right piece exceeds the bank)
        v = np.asarray(v, float); ok = np.isfinite(v)
        return AC.group_means(v[ok], [p for p, o in zip(pdbs, ok) if o]) | {"n": int(ok.sum())}
    agg["m2_mid"] = {k: GN([r["m2"][str(r["m2_mid_cut"])][k] for r in rows]) for k in m2keys}
    agg["m2_bestcut"] = {k: GN([np.nanmin([d[k] for d in r["m2"].values()]) for r in rows]) for k in m2keys}
    agg["m2_meancut"] = {k: GN([np.nanmean([d[k] for d in r["m2"].values()]) for r in rows]) for k in m2keys}
    agg["m3_equal"] = {k: G([r["m3"][r["m3_equal"]][k] for r in rows]) for k in ("argmin_concat", "combo15", "local_err")}
    agg["m3_bestcomp"] = {k: G([min(d[k] for d in r["m3"].values()) for r in rows]) for k in ("argmin_concat", "combo15", "local_err")}
    r4 = [r for r in rows if "m4_equal" in r]
    agg["m4_equal"] = {k: AC.group_means([r["m4_equal"][k] for r in r4], [r["pdb"] for r in r4])
                       for k in ("argmin_concat", "combo8", "local_err", "junction_floor")}
    agg["m4_equal"]["n"] = len(r4)
    # per-target table for the findings
    agg["per_target"] = [{"pdb": r["pdb"], "n": r["n"], "fold": r["fold"], "fail18": r["pdb"] in I.FAIL18,
                          "native_rebuild": r["native_rebuild"], "univ_best_real": r["univ_best_real"],
                          "univ_best_ideal": r["univ_best_ideal"],
                          "m2_mid_argmin": r["m2"][str(r["m2_mid_cut"])]["argmin_concat"],
                          "m2_mid_combo50": r["m2"][str(r["m2_mid_cut"])]["combo50"],
                          "m2_best_combo50": min(d["combo50"] for d in r["m2"].values()),
                          "m2_mid_junction": r["m2_mid_junction_floor"],
                          "m3_equal_combo15": r["m3"][r["m3_equal"]]["combo15"],
                          "m3_best_combo15": min(d["combo15"] for d in r["m3"].values())} for r in rows]
    # paired comparisons vs the whole-window controls
    folds = [r["fold"] for r in rows]
    a = np.array([r["m2"][str(r["m2_mid_cut"])]["combo50"] for r in rows])
    agg["paired_m2mid_combo50_vs_univ_ideal"] = I.paired(a, [r["univ_best_ideal"] for r in rows], folds=folds, names=pdbs)
    agg["paired_m2mid_combo50_vs_native_rebuild"] = I.paired(a, [r["native_rebuild"] for r in rows], folds=folds, names=pdbs)
    b = np.array([r["m3"][r["m3_equal"]]["combo15"] for r in rows])
    agg["paired_m3eq_combo15_vs_m2mid_combo50"] = I.paired(b, a, folds=folds, names=pdbs)
    return agg


if __name__ == "__main__":
    import multiprocessing as mp
    tg = I.targets()
    print("free GB", I.free_gb(), flush=True)
    rows = []
    with mp.Pool(2, initializer=_init) as pool:
        for k, r in enumerate(pool.imap_unordered(run_target, tg)):
            rows.append(r)
            print(f"[{k+1}/126] {r['pdb']} n={r['n']} nat_rebuild={r['native_rebuild']:.3f} univ_ideal={r['univ_best_ideal']:.3f} "
                  f"m2mid_combo={r['m2'][str(r['m2_mid_cut'])]['combo50']:.3f} m3eq_combo={r['m3'][r['m3_equal']]['combo15']:.3f} "
                  f"({r.get('seconds', 0):.0f}s)", flush=True)
    rows.sort(key=lambda r: r["pdb"])
    agg = aggregate(rows)
    for k in range(AC.LMIN, AC.LMAX + 1):
        agg["local"][k] = agg["local"].pop(k)
    print(I.write("assembly_e1_floors", agg))
    for key in ("native_rebuild", "univ_best_real", "univ_best_ideal", "pool_best_real"):
        print(key, agg[key])
    print("m2_mid", agg["m2_mid"]); print("m2_bestcut", agg["m2_bestcut"]); print("m3_equal", agg["m3_equal"])
    print("m3_bestcomp", agg["m3_bestcomp"]); print("m4_equal", agg["m4_equal"])
    print("junction floors", agg["m2_mid_junction_floor"], agg["m3_equal_junction_floor"])
