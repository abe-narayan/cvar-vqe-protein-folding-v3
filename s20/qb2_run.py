"""SPRINT 20 / WORKSTREAM B -- the runner for Q1 (the LEGACY vs AMBER landscape question).

    python -m s20.qb2_run land   [n_targets]     the landscape panel        (PREREG section 3)
    python -m s20.qb2_run opt    [n_targets]     S1, the optimiser battery  (PREREG section 4)
    python -m s20.qb2_run cvar   [n_targets]     S2, the CVaR alpha panel
    python -m s20.qb2_run enc    [n_targets]     S3, theta vs (cos, sin)
    python -m s20.qb2_run all    [n_targets]     all four, in that order

ONLY THE HAMILTONIAN CHANGES.  Sequence, representation, candidate set, ansatz, qubit count,
initialisation, optimiser, evaluation budget, CVaR rule, measurement budget and convergence
threshold are identical across `LEG`, `AMB`, `AMBc` and `DIST` by construction: the same `z0`
array, the same `stable_rng` key and the same `B` are handed to every arm of every objective.

Checkpointed per target; re-running resumes.  Every Ca-RMSD is ORACLE and post-hoc.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

from s20 import qb2_lib as L
from s20 import qb2_opt as OP
from s15 import seed as SD

KINDS = ("LEG", "AMB", "AMBc", "DIST")
B_OPT = 512           # evaluation budget per (target, objective, arm, seed)
SEEDS = (0, 1)
R_START = 7           # starts for the landscape panel (6 basin draws + the cached start)
N_HESS = 2            # Hessians per (target, objective) -- LOCAL, never global topology
LINE_DIRS = 3


def log(tag, msg):
    p = os.path.join(L.RESULTS, f"qb2_{tag}.log")
    with open(p, "a") as fh:
        fh.write(msg + "\n")
    print(msg, flush=True)


def _ck_path(tag):
    return os.path.join(L.RESULTS, f"qb2_{tag}.json")


def ck_load(tag):
    p = _ck_path(tag)
    if os.path.exists(p):
        with open(p) as fh:
            return json.load(fh)
    return {}


def ck_save(tag, d):
    p = _ck_path(tag)
    tmp = p + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(d, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    for _ in range(30):
        try:
            os.replace(tmp, p); return
        except PermissionError:
            time.sleep(0.4)


def mem_hold(min_gb=0.8, tag=""):
    from s12 import instrument as I
    waited = 0
    while I.free_gb() < min_gb and waited < 600:
        log(tag, f"[wait] free RAM {I.free_gb():.2f} GB")
        time.sleep(30); waited += 30


# ============================================================ starts, shared by every arm
def starts_for(tgt):
    """The IDENTICAL start set for every objective and every arm of that target."""
    rng = SD.stable_rng(tgt["pdb"], "starts", salt=L.SALT)
    Z = L.basin_starts(tgt, R_START - 1, rng)
    cs = L.cached_start(tgt["pdb"])
    if cs is not None:
        Z = np.concatenate([cs[None], Z], 0)
    else:
        Z = np.concatenate([L.basin_starts(tgt, 1, rng), Z], 0)
    return Z[:R_START]


def prep(tgt, kinds=KINDS):
    """Build the Hamiltonians, run the BIT-EXACTNESS gate, and do the FREE standardisation."""
    H, sp = L.hams(tgt, kinds, budget=10 ** 9, keep=False)
    cal = {}
    ver = None
    if sp is not None:
        ver = sp.verify(tgt["PHI"], tgt["PSI"], m=4)
    # AMBER on the K=500 pool is computed ONCE and shared by AMB and AMBc (AMBc is a
    # monotone function of the same numbers), instead of paying 500 single points twice.
    e_amb = sp.batch(tgt["PHI"], tgt["PSI"]) if sp is not None else None
    if sp is not None:
        sp.pool_e = e_amb          # cached so the ranking axis costs no extra AMBER calls
    for k in kinds:
        pre = None
        if k == "AMB" and e_amb is not None:
            pre = e_amb
        elif k == "AMBc" and e_amb is not None:
            pre = np.sign(e_amb) * np.log1p(np.abs(e_amb))
        cal[k] = H[k].calibrate(tgt["PHI"], tgt["PSI"], e=pre)
    return H, sp, cal, ver


# ============================================================ MODE: land
def run_land(pdbs):
    tag = "land"
    out = ck_load(tag)
    for ti, pdb in enumerate(pdbs):
        if pdb in out:
            continue
        mem_hold(tag=tag)
        t0 = time.time()
        tgt = L.target(pdb)
        n = tgt["n"]
        Z = starts_for(tgt)
        H, sp, cal, ver = prep(tgt)
        rec = {"n": n, "fold": tgt["fold"], "cal": cal,
               "amber_exact_maxrel": ver[0] if ver else None,
               "amber_exact_ncmp": ver[1] if ver else None,
               "starts_rmsd_ORACLE": L.rmsd_of(Z, tgt).tolist(), "obj": {}}
        for k in KINDS:
            h = H[k]
            rng = SD.stable_rng(pdb, "land", k, salt=L.SALT)
            # the gradient is of the STANDARDISED field: grad(E)/sd_ref.
            gs = []
            for r in range(R_START):
                gr = L.fd_grad(h, Z[r], n, budgeted=False)
                gs.append(float(np.linalg.norm(gr) / h.iqr_ref)
                          if gr is not None and np.isfinite(gr).all() else float("nan"))
            spec = [L.spectrum(L.fd_hess(h, Z[r], n), h.iqr_ref) for r in range(N_HESS)]
            lines = []
            for r in range(2):
                for _d in range(LINE_DIRS):
                    ls = L.line_scan(h, Z[r], n, rng)
                    if ls:
                        lines.append(ls)
            bars = [L.barrier(h, Z[a], Z[b], n) for a, b in ((0, 1), (1, 2), (2, 3))]
            e_starts = h.raw(*L.unpack(Z, n))
            rec["obj"][k] = {
                "gnorm_std": gs,
                "gnorm_std_mean": float(np.nanmean(gs)),
                "gnorm_std_sd": float(np.nanstd(gs)),
                "spec": spec,
                "lines": lines,
                "barrier": [float(b) for b in bars],
                "E_at_starts": [float(x) for x in e_starts],
                "n_nonfinite_starts": int((~np.isfinite(e_starts)).sum()),
            }
        out[pdb] = rec
        ck_save(tag, out)
        log(tag, f"[{ti+1}/{len(pdbs)}] {pdb} n={n} {time.time()-t0:.0f}s "
                 + " ".join(f"{k}:g={rec['obj'][k]['gnorm_std_mean']:.3g}" for k in KINDS))
    with open(os.path.join(L.RESULTS, "qb2_land_COMPLETE"), "w") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S"))
    return out


# ============================================================ MODE: opt
def _run_arm(kind, tgt, arm, seed, z0, budget, sp, cal_k, vqe=None, embed=False):
    """One budgeted arm.  Returns the pre-registered outcome row."""
    h = L.Ham(kind, tgt, budget=budget, keep=False, sp=sp)
    h.mu_ref = cal_k["mu_ref"]; h.sd_ref = cal_k["sd_ref"]
    h.med_ref = cal_k["med_ref"]; h.iqr_ref = cal_k["iqr_ref"]
    F = OP.Field(h, tgt)
    # the sentinel penalty, identical for every arm of this objective
    F.set_pen(cal_k.get("ref_std_max", 6.0))
    rng = SD.stable_rng(tgt["pdb"], kind, arm, seed, salt=L.SALT)
    t0 = time.time()
    if embed:
        info = _run_embedded(F, z0, rng, arm, vqe)
    elif vqe is not None:
        info = OP.arm_vqe(F, z0, rng, **vqe)
    else:
        info = OP.ARMS[arm](F, z0, rng)
    rm, braw = float("nan"), float("nan")
    if F.best_z is not None:
        # `Field.best_z` is always in TORSION space -- the embedded field delegates its
        # __call__ to the inner Field after the arctan2 retraction, so no re-mapping is needed.
        rm = float(L.rmsd_of(F.best_z[None], tgt)[0])
        braw = float(h.raw(*L.unpack(F.best_z[None], tgt["n"]))[0])
    return {"best_std": float(F.best), "best_raw": braw, "rmsd_ORACLE": rm, "used": int(h.used),
            "n_nonfinite": int(F.n_nonfinite), "wall": time.time() - t0, **info}


class _EmbField:
    """S3: the SAME objective seen through u = (cos theta, sin theta), with a unit-norm
    retraction.  Physically identical, geometrically different.  No discretisation."""

    def __init__(self, F):
        self.F = F
        self.nz = int(F.n)          # residues
        self.n = 2 * int(F.n)       # HALF the embedded dimension: the arms size FD as 2*n
        self.t = F.t

    @staticmethod
    def wrap(u):
        return u

    @property
    def left(self):
        return self.F.left

    def restart(self, rng, k=1):
        Z = self.F.restart(rng, k)
        return np.concatenate([np.cos(Z), np.sin(Z)], axis=1)

    def to_z(self, U):
        """The retraction.  `arctan2` is scale-invariant, so the radius is a pure gauge: the
        embedded coordinate is physically equivalent to theta and geometrically different."""
        U = np.atleast_2d(np.asarray(U, float))
        d = 2 * self.nz
        return np.arctan2(U[:, d:], U[:, :d])

    def __call__(self, U):
        return self.F(self.to_z(U))


def _run_embedded(F, z0, rng, arm, vqe):
    E = _EmbField(F)
    u0 = np.concatenate([np.cos(z0), np.sin(z0)])
    return OP.ARMS[arm](E, u0, rng) if arm in OP.ARMS else {}


def run_opt(pdbs, arms=None, tag="opt", budget=B_OPT, kinds=KINDS, embed=False):
    arms = arms or list(OP.ARMS)
    out = ck_load(tag)
    for ti, pdb in enumerate(pdbs):
        if pdb in out:
            continue
        mem_hold(tag=tag)
        t0 = time.time()
        tgt = L.target(pdb)
        Z = starts_for(tgt)
        H, sp, cal, ver = prep(tgt, kinds)
        rec = {"n": tgt["n"], "fold": tgt["fold"], "cal": cal,
               "amber_exact_maxrel": ver[0] if ver else None,
               "budget": budget, "rows": {}}
        for k in kinds:
            for arm in arms:
                for sd_ in SEEDS:
                    z0 = Z[sd_ % len(Z)]
                    r = _run_arm(k, tgt, arm, sd_, z0, budget, sp, cal[k], embed=embed)
                    rec["rows"][f"{k}|{arm}|{sd_}"] = r
        # zero-information RMSD references, ORACLE, 0 evaluations
        # THE RANKING AXIS, free: does this objective order the SAME 500 candidates the
        # pipeline actually uses?  Connects every landscape metric to structural outcome.
        Zpool = L.pack(tgt["PHI"], tgt["PSI"])
        rr = L.rmsd_of(Zpool, tgt)
        rec["rank"] = {}
        for k in kinds:
            if k == "AMB" and sp is not None:
                ep = np.asarray(sp.pool_e, float)
            elif k == "AMBc" and sp is not None:
                ea = np.asarray(sp.pool_e, float)
                ep = np.sign(ea) * np.log1p(np.abs(ea))
            else:
                ep = H[k].raw(tgt["PHI"], tgt["PSI"])
            fin = np.isfinite(ep)
            rec["rank"][k] = {
                "spearman_E_vs_rmsd_ORACLE": float(L.spearman(ep[fin], rr[fin])),
                "argmin_rmsd_ORACLE": float(rr[fin][int(np.argmin(ep[fin]))]),
                "top1pct_mean_rmsd_ORACLE": float(
                    rr[fin][np.argsort(ep[fin])[:max(1, int(0.01 * fin.sum()))]].mean()),
            }
        rec["ref"] = {
            "start_rmsd_ORACLE": float(L.rmsd_of(Z[0][None], tgt)[0]),
            "pool_best_ORACLE": float(np.min(L.rmsd_of(L.pack(tgt["PHI"], tgt["PSI"]), tgt))),
            "pool_mean_ORACLE": float(np.mean(L.rmsd_of(L.pack(tgt["PHI"], tgt["PSI"]), tgt))),
        }
        out[pdb] = rec
        ck_save(tag, out)
        s = " ".join(f"{k}:{np.nanmin([rec['rows'][x]['rmsd_ORACLE'] for x in rec['rows'] if x.startswith(k+'|')]):.2f}"
                     for k in kinds)
        log(tag, f"[{ti+1}/{len(pdbs)}] {pdb} n={tgt['n']} {time.time()-t0:.0f}s bestRMSD {s}")
    with open(os.path.join(L.RESULTS, f"qb2_{tag}_COMPLETE"), "w") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S"))
    return out


# ============================================================ MODE: cvar
ALPHAS = (0.05, 0.25, 1.00)


def run_cvar(pdbs, budget=B_OPT, kinds=("LEG", "AMB", "AMBc", "DIST")):
    tag = "cvar"
    out = ck_load(tag)
    for ti, pdb in enumerate(pdbs):
        if pdb in out:
            continue
        mem_hold(tag=tag)
        t0 = time.time()
        tgt = L.target(pdb)
        Z = starts_for(tgt)
        H, sp, cal, ver = prep(tgt, kinds)
        rec = {"n": tgt["n"], "fold": tgt["fold"], "budget": budget, "rows": {}}
        for k in kinds:
            for a in ALPHAS:
                for sd_ in SEEDS:
                    r = _run_arm(k, tgt, f"vqe{a}", sd_, Z[sd_ % len(Z)], budget, sp, cal[k],
                                 vqe={"alpha": a, "shots": 64, "train": True})
                    rec["rows"][f"{k}|vqe{a}|{sd_}"] = r
            # MANDATORY: the same circuit at theta_0, zero gradient steps, same budget
            for sd_ in SEEDS:
                r = _run_arm(k, tgt, "vqe_untrained", sd_, Z[sd_ % len(Z)], budget, sp, cal[k],
                             vqe={"alpha": 1.0, "shots": 64, "train": False})
                rec["rows"][f"{k}|vqe_untrained|{sd_}"] = r
        out[pdb] = rec
        ck_save(tag, out)
        log(tag, f"[{ti+1}/{len(pdbs)}] {pdb} n={tgt['n']} {time.time()-t0:.0f}s")
    with open(os.path.join(L.RESULTS, "qb2_cvar_COMPLETE"), "w") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S"))
    return out


# ============================================================ main
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "land"
    n_t = int(sys.argv[2]) if len(sys.argv) > 2 else L.N_SUBSET
    ss = L.subset(L.N_SUBSET)[:n_t]
    pdbs = [t["pdb"] for t in ss]
    L.write("qb2_config", {"salt": L.SALT, "subset": pdbs, "B_OPT": B_OPT, "SEEDS": list(SEEDS),
                           "R_START": R_START, "FD_H": L.FD_H, "N_HESS": N_HESS,
                           "ALPHAS": list(ALPHAS), "kinds": list(KINDS),
                           "mode": mode, "n_targets": len(pdbs)}, complete=True)
    # S2 and S3 run on a declared PREFIX of the same subset (compute-bound); every claim
    # carries its own n and the prefix is persisted in `qb2_config.json`.
    n_small = min(len(pdbs), 10)
    if mode in ("land", "all"):
        run_land(pdbs)
    if mode in ("opt", "all"):
        run_opt(pdbs)
    if mode in ("cvar", "all"):
        run_cvar(pdbs[:n_small])
    if mode in ("enc", "all"):
        run_opt(pdbs[:n_small], arms=["spsa", "adam_fd", "lbfgs_fd", "nelder"], tag="enc",
                kinds=("LEG", "AMB", "AMBc"), embed=True)
    print("done")


if __name__ == "__main__":
    main()
