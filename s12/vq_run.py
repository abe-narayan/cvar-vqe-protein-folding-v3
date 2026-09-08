"""s12 QUANTUM-ROLE -- stage 1: is the problem genuinely interacting, and how hard is it?

For every VQ-A instance:
  * an EXACT functional-ANOVA decomposition of the energy variance under the uniform
    product measure into one-body and pure-interaction parts (the quantitative form of
    "the value of a choice depends on the other choices");
  * the certified exhaustive optimum;
  * a SEPARABLE control (minimise sum_s h_s only, i.e. pretend the couplings are absent)
    -- how often it finds the true optimum, and what it costs;
  * the number of distinct 1-change local minima (landscape multimodality);
  * every classical solver, with wall time and objective-evaluation count;
  * the emitted CA-RMSD of each solver's answer (evaluation only -- the native is a label).

Usage: python -m s12.vq_run stage1 <family>
"""
from __future__ import annotations
import os, sys, json, time, math
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import vq_lib as V
from s12 import vq_classical as C
from s12.vq_build import FAMILIES

NATCACHE = os.path.join(ROOT, "s12", "cache", "vq", "nat_ca.npz")


def native(pdb):
    """`nat_ca` for a target (ORACLE -- evaluation only), cached compactly."""
    if not hasattr(native, "_z"):
        native._z = dict(np.load(NATCACHE)) if os.path.exists(NATCACHE) else {}
    if pdb not in native._z:
        native._z[pdb] = I.load_univ(pdb)["nat_ca"]
        np.savez_compressed(NATCACHE, **native._z)
    return native._z[pdb]


# ------------------------------------------------------------------ structure measures
def anova(inst):
    """Exact variance decomposition of F under the uniform product measure on choices.

    F(x) = c + sum_s a_s(x_s) + sum_{s<t} b_st(x_s,x_t) with every term centred:
      a_s   = (h_s - mean h_s)/npairs + sum_{t!=s} (row/col means of J_st, centred)
      b_st  = J_st - rowmean - colmean + grandmean, over npairs
    Under independent uniform x_s the terms are orthogonal, so
      Var(F) = sum_s Var(a_s) + sum_{s<t} Var(b_st)     -- checked against a direct
    enumeration of Var(F) in `interaction_share`.
    """
    k, m, npz = inst["k"], inst["m"], inst["npairs"]
    a = [np.asarray(inst["h"][s], float) / npz for s in range(k)]
    a = [x - x.mean() for x in a]
    binter = {}
    for (s, u), Jm in inst["J"].items():
        Jm = np.asarray(Jm, float) / npz
        rm, cm, gm = Jm.mean(1), Jm.mean(0), Jm.mean()
        a[s] = a[s] + (rm - gm)
        a[u] = a[u] + (cm - gm)
        binter[(s, u)] = Jm - rm[:, None] - cm[None, :] + gm
    v1 = float(sum(x.var() for x in a))
    v2 = float(sum(B.var() for B in binter.values()))
    return v1, v2


def interaction_share(inst, E=None):
    if E is None:
        _, E = V.enumerate_all(inst)
    v1, v2 = anova(inst)
    tot = float(np.var(E))
    return dict(var_total=tot, var_onebody=v1, var_interaction=v2,
                anova_residual=float(abs(tot - v1 - v2)),
                interaction_share=float(v2 / max(tot, 1e-300)))


def separable_opt(inst):
    """argmin of the one-body part alone (the couplings pretended absent)."""
    return np.array([int(np.argmin(inst["h"][s])) for s in range(inst["k"])], int)


def separable_anova(inst):
    """argmin of the BEST ADDITIVE approximation of F (the ANOVA one-body terms).

    This is the strongest possible "ignore the interactions" control: it is not merely the
    one-body table h, it is the optimum of the closest separable function to F in L2 under
    the uniform product measure.  If the optimum of THIS still misses the true optimum, the
    coupling is doing real work.
    """
    k, m, npz = inst["k"], inst["m"], inst["npairs"]
    a = [np.asarray(inst["h"][s], float) / npz for s in range(k)]
    for (s, u), Jm in inst["J"].items():
        Jm = np.asarray(Jm, float) / npz
        a[s] = a[s] + Jm.mean(1)
        a[u] = a[u] + Jm.mean(0)
    return np.array([int(np.argmin(x)) for x in a], int)


def n_local_minima(cfgs, E, k, m):
    """Number of configurations that are minima of the 1-change neighbourhood."""
    Eg = E.reshape([m] * k)
    is_min = np.ones(Eg.shape, bool)
    for s in range(k):
        mn = Eg.min(axis=s, keepdims=True)
        is_min &= (Eg <= mn + 1e-12)
    return int(is_min.sum())


# ------------------------------------------------------------------ stage 1
def stage1(fam, limit=None, do_milp=True, out=None):
    cfg = FAMILIES[fam]
    k, m = cfg["k"], cfg["m"]
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    rows = []
    t00 = time.time()
    for q, t in enumerate(tg):
        if V.segments(t["n"], k, cfg["minlen"]) is None:
            continue
        inst = V.cached_instance(t, **cfg)
        h, J, npz = inst["h"], inst["J"], inst["npairs"]
        nat = native(t["pdb"])

        exh = C.a_exhaustive(h, J, npz)
        E, cfgs = exh.pop("all_E"), exh.pop("all_cfg")
        opt = exh["value"]
        st = interaction_share(inst, E)
        st["n_local_minima"] = n_local_minima(cfgs, E, k, m)
        st["n_configs"] = int(len(E))

        res = {"exhaustive": exh}
        sep = separable_opt(inst)
        res["separable_h"] = dict(cfg=sep.tolist(), value=float(C.a_energy(h, J, sep, npz)[0]),
                                  evals=k * m, wall=0.0)
        sepa = separable_anova(inst)
        res["separable_anova"] = dict(cfg=sepa.tolist(),
                                      value=float(C.a_energy(h, J, sepa, npz)[0]),
                                      evals=k * m, wall=0.0)
        res["greedy"] = C.a_greedy(h, J, npz)
        res["greedy_ls"] = C.a_greedy_ls(h, J, npz)
        res["beam4"] = C.a_beam(h, J, npz, width=4)
        res["beam16"] = C.a_beam(h, J, npz, width=16)
        res["anneal1k"] = C.a_anneal(h, J, npz, steps=1000, seed=0)
        res["anneal_matched"] = C.a_anneal(h, J, npz, steps=max(64, res["greedy_ls"]["evals"]), seed=0)
        res["random_matched"] = C.a_random(h, J, npz, n=res["greedy_ls"]["evals"], seed=0)
        if do_milp:
            try:
                res["milp"] = C.a_milp(h, J, npz)
            except Exception as e:                     # noqa: BLE001
                res["milp"] = dict(error=str(e))

        # accuracy of every answer (evaluation only)
        rm = {}
        for nm, r in res.items():
            if "cfg" not in r:
                continue
            rm[nm] = float(I.ca_rmsd(V.assemble(inst, r["cfg"]), nat))
        # oracle reference points inside this instance family
        Xall = V.assemble_batch(inst, cfgs)
        rr = I.kabsch_rmsd_batch(Xall, nat)
        rm["ORACLE_best_assembly"] = float(rr.min())
        rm["mean_assembly"] = float(rr.mean())
        rm["anchor_fit_ca"] = float(I.ca_rmsd(inst["anchor"], nat))
        rank_of_opt = float((rr < rr[int(np.argmin(E))]).mean())
        rows.append(dict(pdb=t["pdb"], n=t["n"], fold=t["fold"], fam=fam,
                         structure=st, solvers=res, rmsd=rm,
                         rmsd_percentile_of_optimum=rank_of_opt,
                         corr_E_rmsd=float(np.corrcoef(E, rr)[0, 1])))
        if q % 20 == 0:
            print(f"  [{fam}] {q+1}/{len(tg)} {time.time()-t00:.0f}s free={I.free_gb():.2f}", flush=True)
    path = I.write(out or f"vq_stage1_{fam}", dict(family=fam, cfg=cfg, rows=rows))
    print("wrote", path, flush=True)
    return rows


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "stage1":
        fam = sys.argv[2]
        lim = int(sys.argv[3]) if len(sys.argv) > 3 else None
        stage1(fam, lim, do_milp=(FAMILIES[fam]["m"] <= 64))
