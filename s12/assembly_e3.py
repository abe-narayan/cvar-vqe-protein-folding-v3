"""E3 -- DEPLOYABLE assembly as a combinatorial problem, and the assembled pool through the
shipped terminal operator.

Per target (deployable view only, natives used for evaluation):
  1. cut structures: every 2-piece cut with both lengths in [4,10]; every 3-piece composition
     with lengths in [4,8].
  2. at every piece position the top-k pieces by the E2 key (default `combo`), pieces whose
     junction-relevant torsions are chain-end placeholders excluded.
  3. every combination assembled by torsion concatenation (cut convention) and built.
  4. objectives (lower = better):
        dist  = shipped distogram Bayes-risk score of the CA trace
        junc  = sum over junctions of -log P(state_{b+1} | state_b)  (library bigram, per fold)
        rg    = |rg - rg_expected(n)| / 1.0   (rg_expected from the fold's library windows)
        legacy= Legacy total energy of the ideal backbone (top-2000 by dist only)
     combined objectives are sums of per-target z-scores.
  5. solvers: exhaustive (ground truth), beam (left-to-right, width B), simulated annealing.
  6. outputs: argmin, top-M candidate pool (M=500), the pool through the terminal operator
     (distogram top-75 -> coordinate_average -> project), hybrid with the retrieval pool.
"""
from __future__ import annotations
import os, sys, json, time, math
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import assembly_bank as AB
from s12 import assembly_common as AC
from s12 import assembly_e2 as E2

KEY = os.environ.get("ASM_KEY", "combo")
KTOP = int(os.environ.get("ASM_K", "20"))
M_POOL = 500
DO_PROJECT = os.environ.get("ASM_PROJECT", "1") == "1"
L2 = (4, 10); L3 = (4, 8)


def structures(n):
    out = [(c, n - c) for c in range(L2[0], n - L2[0] + 1) if c <= L2[1] and n - c <= L2[1]]
    lo = L3[0] if n >= 3 * L3[0] else 3      # n < 12 has no 3-piece split with L >= 4
    out += AC.cuts_m(n, 3, lo, L3[1])
    return out


def valid_junction_mask(fold, L, first, last):
    R = AB.residues(fold); P = AB.pieces(fold, L); st = P["start"]
    ok = np.ones(len(st), bool)
    if not first:
        ok &= R["valid_phi"][st]
    if not last:
        ok &= R["valid_psi"][st + L - 1]
    return ok


_TOP = {}


def top_pieces(t, K, s, L, first, last, k):
    key = (t["pdb"], KEY, s, L, first, last, k)
    if key not in _TOP:
        sc = E2.piece_scores(t, K, s, L, only=KEY)[KEY].astype(float)
        sc[~valid_junction_mask(t["fold"], L, first, last)] = -np.inf
        _TOP[key] = np.argsort(-sc, kind="stable")[:k]
    return _TOP[key]


def enumerate_structure(t, K, comp, k):
    """All k^m assemblies for one composition: returns PHI, PSI (N, n), choice (N, m), piece ids."""
    n, fold = t["n"], t["fold"]; iv = AC.intervals(comp); m = len(iv)
    tops, tors, states = [], [], []
    R = AB.residues(fold)
    for q, (s, L) in enumerate(iv):
        top = top_pieces(t, K, s, L, q == 0, q == m - 1, k)
        ph, ps, _ = AB.piece_torsions(fold, L, top)
        idx = AB.pieces(fold, L)["start"][top]
        tops.append(top); tors.append((ph, ps)); states.append((R["state"][idx], R["state"][idx + L - 1]))
    grids = np.array(np.meshgrid(*[np.arange(len(tp)) for tp in tops], indexing="ij")).reshape(m, -1).T
    N = len(grids)
    PHI = np.zeros((N, n)); PSI = np.zeros((N, n))
    for q, (s, L) in enumerate(iv):
        PHI[:, s:s + L] = tors[q][0][grids[:, q]]; PSI[:, s:s + L] = tors[q][1][grids[:, q]]
    junc = np.zeros(N)
    T = R["trans"]
    for q in range(m - 1):
        a = states[q][1][grids[:, q]]; b = states[q + 1][0][grids[:, q + 1]]
        junc += -np.log(T[a, b])
    return PHI, PSI, grids, tops, junc


def rg_expected(u):
    return float(np.median(AC.rg(u["W"][I.pool_idx(u)])))


def legacy_total(seq, PHI, PSI):
    from core import geometry as geo
    from core import energy as en
    bb = geo.build_backbone_batch(np.asarray(PHI, float), np.asarray(PSI, float))
    comp = en.components_batch(seq, bb, phi=np.asarray(PHI, float), psi=np.asarray(PSI, float))
    return sum(en.DEFAULT_WEIGHTS[k] * comp[k] for k in en.TERM_NAMES)


def z(v):
    v = np.asarray(v, float); return (v - v.mean()) / (v.std() + 1e-9)


def beam_search(t, K, comp, k, width, dg, i, j):
    """Left-to-right beam over positions; partial assemblies scored by the distogram score of the
    pairs already built (a deployable, monotone-ish surrogate). Returns the beam's best full choice."""
    n, fold = t["n"], t["fold"]; iv = AC.intervals(comp); m = len(iv)
    tops, tors = [], []
    for q, (s, L) in enumerate(iv):
        top = top_pieces(t, K, s, L, q == 0, q == m - 1, k)
        ph, ps, _ = AB.piece_torsions(fold, L, top); tops.append(top); tors.append((ph, ps))
    beam = [((), np.zeros(0), np.zeros(0))]
    for q, (s, L) in enumerate(iv):
        cand = []
        PHI = []; PSI = []; parents = []
        for (ch, ph, ps) in beam:
            for a in range(len(tops[q])):
                PHI.append(np.concatenate([ph, tors[q][0][a]])); PSI.append(np.concatenate([ps, tors[q][1][a]])); parents.append((ch + (a,)))
        PHI = np.array(PHI); PSI = np.array(PSI)
        end = s + L
        ca = AC.build_many(PHI, PSI)
        sel = (i < end) & (j < end)
        if sel.sum() > 0:
            D = I.pair_dists(ca, i[sel], j[sel])
            grid = dg["grid"]; risk = dg["risk"][sel]
            g = np.clip(((D - grid[0]) / 0.05).astype(np.int32), 0, len(grid) - 1)
            sc = risk[np.arange(risk.shape[0])[None, :], g].mean(1)
        else:
            sc = np.zeros(len(PHI))
        order = np.argsort(sc, kind="stable")[:width]
        beam = [(parents[o], PHI[o], PSI[o]) for o in order]
    return beam[0][0]


def anneal(obj_lookup, shape, steps=4000, seed=0):
    """Simulated annealing over the choice vector; obj_lookup(choice tuple) -> energy (from the
    exhaustive table, so this measures the SOLVER, not the objective)."""
    rng = np.random.default_rng(seed); m = len(shape)
    x = tuple(int(rng.integers(0, shape[q])) for q in range(m)); e = obj_lookup(x); best = (e, x)
    for step in range(steps):
        T = 1.0 * (1 - step / steps) + 0.02
        q = int(rng.integers(0, m)); y = list(x); y[q] = int(rng.integers(0, shape[q])); y = tuple(y)
        e2 = obj_lookup(y)
        if e2 < e or rng.random() < math.exp(-(e2 - e) / T):
            x, e = y, e2
            if e < best[0]:
                best = (e, x)
    return best


def run_target(t):
    pdb, n, fold, seq = t["pdb"], t["n"], t["fold"], t["seq"]
    out_path = os.path.join(I.CACHE, f"asm_e3_{KEY}_k{KTOP}_{pdb}.json")
    if os.path.exists(out_path):
        return json.load(open(out_path))
    t0 = time.time()
    u = I.load_univ(pdb); nat = u["nat_ca"]; p = I.pool_idx(u)
    dg = I.distogram(pdb); i, j = I.pair_index(n)
    K = E2.target_keys(t); tk = time.time()
    # ---- enumerate every structure
    allPHI, allPSI, allJ, tag, meta = [], [], [], [], []
    for comp in structures(n):
        PHI, PSI, grids, tops, junc = enumerate_structure(t, K, comp, KTOP)
        allPHI.append(PHI); allPSI.append(PSI); allJ.append(junc); tag += [len(meta)] * len(PHI)
        meta.append({"comp": comp, "grids": grids, "tops": tops})
    PHI = np.concatenate(allPHI); PSI = np.concatenate(allPSI); junc = np.concatenate(allJ); tag = np.array(tag)
    N = len(PHI); te = time.time()
    CA = AC.build_many(PHI, PSI); tb = time.time()
    D = I.pair_dists(CA, i, j)
    dist = I.shipped_score(dg, D)
    rr = I.kabsch_rmsd_batch(CA, nat)                 # ORACLE, evaluation only
    rgs = AC.rg(CA); rge = rg_expected(u); rgpen = np.abs(rgs - rge)
    # legacy on the top-2000 by dist
    top_d = np.argsort(dist, kind="stable")[:2000]
    leg = np.full(N, np.nan); leg[top_d] = legacy_total(seq, PHI[top_d], PSI[top_d])
    objs = {"dist": dist, "dist+junc": z(dist) + z(junc), "dist+junc+rg": z(dist) + z(junc) + z(rgpen)}
    lz = np.full(N, np.inf); lz[top_d] = z(dist[top_d]) + z(junc[top_d]) + z(leg[top_d])
    objs["dist+junc+legacy"] = lz
    res = {"pdb": pdb, "n": n, "fold": fold, "n_candidates": int(N), "n_structures": len(meta),
           "gen_best_rr": float(rr.min()), "gen_median_rr": float(np.median(rr)),
           "gen_frac_under2": float((rr < 2.0).mean()), "gen_frac_under1": float((rr < 1.0).mean()),
           "pool_best_rr": float(u["rr"][p].min()), "univ_best_rr": float(u["rr"].min()),
           "rg_expected": rge, "rg_native": float(AC.rg(nat)), "arms": {}}
    # per-structure best (which cut structure generates best) -- oracle diagnostic
    res["structure_best_rr"] = {"-".join(map(str, meta[s]["comp"])): float(rr[tag == s].min()) for s in range(len(meta))}
    mm = lambda c: float(min([v for kk, v in res["structure_best_rr"].items() if kk.count("-") == c], default=float("nan")))
    res["m2_gen_best_rr"] = mm(1); res["m3_gen_best_rr"] = mm(2)
    # retrieval pool objects for hybrid / control
    Wp = u["W"][p]; Dp = I.pair_dists(Wp, i, j); sp = I.shipped_score(dg, Dp.astype(np.float32).astype(float)); rp = u["rr"][p]
    rec = I.shipped_record(pdb)
    res["shipped_fit_rr"] = I.ca_rmsd(np.asarray(rec["fit_ca"]), nat)
    res["shipped_ca_rr"] = I.ca_rmsd(np.asarray(rec["ca"]), nat)
    res["shipped_argmin_rr"] = float(rp[int(np.argmin(sp))])
    for name, ob in objs.items():
        order = np.argsort(ob, kind="stable")
        am = int(order[0]); top = order[:M_POOL]
        arm = {"argmin_rr": float(rr[am]), "argmin_dist": float(dist[am]), "pool_best_rr": float(rr[top].min()),
               "pool_mean_rr": float(rr[top].mean()), "argmin_structure": "-".join(map(str, meta[tag[am]]["comp"]))}
        # terminal operator on the assembled pool: distogram top-75 -> avg -> project
        sub = top[np.argsort(dist[top], kind="stable")[:I.M]]
        arm["top75_best_rr"] = float(rr[sub].min())
        C, _ = I.coordinate_average(CA[sub]); arm["avg_rr"] = I.ca_rmsd(C, nat)
        # hybrid: retrieval 500 U assembled 500
        Wh = np.concatenate([Wp, CA[top]]); sh = np.concatenate([sp, dist[top]]); rh = np.concatenate([rp, rr[top]])
        subh = np.argsort(sh, kind="stable")[:I.M]
        arm["hybrid_pool_best_rr"] = float(rh.min()); arm["hybrid_top75_best_rr"] = float(rh[subh].min())
        arm["hybrid_top75_frac_assembled"] = float((subh >= len(Wp)).mean())
        arm["hybrid_argmin_rr"] = float(rh[int(np.argmin(sh))])
        Ch, _ = I.coordinate_average(Wh[subh]); arm["hybrid_avg_rr"] = I.ca_rmsd(Ch, nat)
        if DO_PROJECT and name in ("dist",):
            pr = I.project(C, seq, fold); arm["proj_fit_rr"] = I.ca_rmsd(pr["fit_ca"], nat); arm["proj_ca_rr"] = I.ca_rmsd(pr["ca"], nat)
            prh = I.project(Ch, seq, fold); arm["hybrid_proj_fit_rr"] = I.ca_rmsd(prh["fit_ca"], nat); arm["hybrid_proj_ca_rr"] = I.ca_rmsd(prh["ca"], nat)
        res["arms"][name] = arm
    # ---- solver validation on the dist+junc objective, per structure: exhaustive vs beam vs SA
    ob = objs["dist+junc"]; solv = []
    for sidx, mt in enumerate(meta):
        rows = np.where(tag == sidx)[0]; sub_ob = ob[rows]
        ex = int(np.argmin(sub_ob)); ex_val = float(sub_ob[ex])
        grids = mt["grids"]; shape = [len(tp) for tp in mt["tops"]]
        lut = {tuple(g): float(v) for g, v in zip(grids, sub_ob)}
        sa_val, _ = anneal(lambda x: lut[x], shape, steps=2000)
        bm = beam_search(t, K, mt["comp"], KTOP, 20, dg, i, j)
        bm_val = lut.get(tuple(bm), float("nan"))
        solv.append({"comp": "-".join(map(str, mt["comp"])), "size": int(len(rows)), "exhaustive": ex_val, "sa2000": sa_val, "beam20": bm_val,
                     "sa_gap_rank": int((sub_ob < sa_val).sum()), "beam_gap_rank": int((sub_ob < bm_val).sum()) if np.isfinite(bm_val) else None})
    res["solvers"] = solv
    res["seconds"] = time.time() - t0
    res["timing"] = {"keys": tk - t0, "enumerate": te - tk, "build": tb - te, "rest": time.time() - tb}
    json.dump(res, open(out_path, "w"), default=float)
    return json.load(open(out_path))


def _init():
    os.environ["OMP_NUM_THREADS"] = "2"
    try:
        import torch; torch.set_num_threads(2)
    except Exception:
        pass


def aggregate(rows):
    pdbs = [r["pdb"] for r in rows]; folds = [r["fold"] for r in rows]
    G = lambda v: AC.group_means(v, pdbs)
    agg = {"key": KEY, "k": KTOP, "n_targets": len(rows), "n_candidates_mean": float(np.mean([r["n_candidates"] for r in rows]))}
    for key in ("gen_best_rr", "gen_median_rr", "gen_frac_under2", "gen_frac_under1", "m2_gen_best_rr", "m3_gen_best_rr",
                "pool_best_rr", "univ_best_rr", "shipped_fit_rr", "shipped_ca_rr", "shipped_argmin_rr"):
        v = np.array([r[key] for r in rows], float); ok = np.isfinite(v)
        agg[key] = AC.group_means(v[ok], [p for p, o in zip(pdbs, ok) if o])
    agg["arms"] = {}
    for name in rows[0]["arms"]:
        agg["arms"][name] = {k: G([r["arms"][name][k] for r in rows]) for k in rows[0]["arms"][name] if not isinstance(rows[0]["arms"][name][k], str)}
    P = {}
    P["asm_pool_best_vs_retrieval_pool_best"] = I.paired([r["arms"]["dist"]["pool_best_rr"] for r in rows], [r["pool_best_rr"] for r in rows], folds=folds, names=pdbs)
    P["asm_gen_best_vs_univ_best"] = I.paired([r["gen_best_rr"] for r in rows], [r["univ_best_rr"] for r in rows], folds=folds, names=pdbs)
    P["asm_argmin_vs_shipped_argmin"] = I.paired([r["arms"]["dist"]["argmin_rr"] for r in rows], [r["shipped_argmin_rr"] for r in rows], folds=folds, names=pdbs)
    P["asm_avg_vs_shipped_fit"] = I.paired([r["arms"]["dist"]["avg_rr"] for r in rows], [r["shipped_fit_rr"] for r in rows], folds=folds, names=pdbs)
    P["hybrid_avg_vs_shipped_fit"] = I.paired([r["arms"]["dist"]["hybrid_avg_rr"] for r in rows], [r["shipped_fit_rr"] for r in rows], folds=folds, names=pdbs)
    P["hybrid_top75_best_vs_retrieval_top75_best(2.306)"] = I.paired([r["arms"]["dist"]["hybrid_top75_best_rr"] for r in rows], [I.shipped_record(r["pdb"]) and None or 0 for r in rows], folds=folds, names=pdbs) if False else None
    if "proj_fit_rr" in rows[0]["arms"]["dist"]:
        for name in ("dist",):
            P[f"{name}_proj_fit_vs_shipped_fit"] = I.paired([r["arms"][name]["proj_fit_rr"] for r in rows], [r["shipped_fit_rr"] for r in rows], folds=folds, names=pdbs)
            P[f"{name}_hybrid_proj_fit_vs_shipped_fit"] = I.paired([r["arms"][name]["hybrid_proj_fit_rr"] for r in rows], [r["shipped_fit_rr"] for r in rows], folds=folds, names=pdbs)
            P[f"{name}_proj_ca_vs_shipped_ca"] = I.paired([r["arms"][name]["proj_ca_rr"] for r in rows], [r["shipped_ca_rr"] for r in rows], folds=folds, names=pdbs)
        a = np.array([r["arms"]["dist"]["proj_fit_rr"] for r in rows]); b = np.array([r["shipped_fit_rr"] for r in rows])
        agg["error_corr_asm_vs_shipped"] = float(np.corrcoef(a, b)[0, 1])
        agg["error_corr_hybrid_vs_shipped"] = float(np.corrcoef([r["arms"]["dist"]["hybrid_proj_fit_rr"] for r in rows], b)[0, 1])
        agg["oracle_min_of_two"] = G(np.minimum(a, b))
    agg["paired"] = {k: v for k, v in P.items() if v is not None}
    solv = [s for r in rows for s in r["solvers"]]
    agg["solvers"] = {"n_problems": len(solv), "mean_size": float(np.mean([s["size"] for s in solv])),
                      "sa_exact_frac": float(np.mean([s["sa2000"] <= s["exhaustive"] + 1e-9 for s in solv])),
                      "sa_mean_rank_gap": float(np.mean([s["sa_gap_rank"] for s in solv])),
                      "beam_exact_frac": float(np.mean([s["beam20"] <= s["exhaustive"] + 1e-9 for s in solv if s["beam_gap_rank"] is not None])),
                      "beam_mean_rank_gap": float(np.mean([s["beam_gap_rank"] for s in solv if s["beam_gap_rank"] is not None]))}
    agg["per_target"] = [{"pdb": r["pdb"], "fail18": r["pdb"] in I.FAIL18, "gen_best": r["gen_best_rr"], "asm_pool_best": r["arms"]["dist"]["pool_best_rr"],
                          "retr_pool_best": r["pool_best_rr"], "asm_argmin": r["arms"]["dist"]["argmin_rr"], "shipped_argmin": r["shipped_argmin_rr"],
                          "asm_proj": r["arms"]["dist"].get("proj_fit_rr"), "hybrid_proj": r["arms"]["dist"].get("hybrid_proj_fit_rr"), "shipped_fit": r["shipped_fit_rr"]} for r in rows]
    return agg


if __name__ == "__main__":
    import multiprocessing as mp
    tg = I.targets()
    print(f"key={KEY} k={KTOP} project={DO_PROJECT} free GB {I.free_gb():.2f}", flush=True)
    rows = []
    with mp.Pool(2, initializer=_init) as pool:
        for k, r in enumerate(pool.imap_unordered(run_target, tg)):
            rows.append(r)
            a = r["arms"]["dist"]
            print(f"[{k+1}/126] {r['pdb']} N={r['n_candidates']} gen_best={r['gen_best_rr']:.2f} pool_best={a['pool_best_rr']:.2f} (retr {r['pool_best_rr']:.2f}) "
                  f"argmin={a['argmin_rr']:.2f} avg={a['avg_rr']:.2f} proj={a.get('proj_fit_rr', float('nan')):.2f} hyb={a.get('hybrid_proj_fit_rr', float('nan')):.2f} shipped={r['shipped_fit_rr']:.2f} ({r['seconds']:.0f}s)", flush=True)
    rows.sort(key=lambda r: r["pdb"])
    agg = aggregate(rows)
    print(I.write(f"assembly_e3_{KEY}_k{KTOP}", agg))
    for key in ("gen_best_rr", "m2_gen_best_rr", "m3_gen_best_rr", "pool_best_rr", "univ_best_rr", "shipped_fit_rr", "shipped_argmin_rr"):
        print(key, agg[key])
    for name, d in agg["arms"].items():
        print(name, {k: round(v["all"], 3) for k, v in d.items()})
    print("solvers", agg["solvers"])
    for k, v in agg["paired"].items():
        print(k, round(v["mean_diff"], 3), v["ci95"], v["n_better"], v["n_worse"])
