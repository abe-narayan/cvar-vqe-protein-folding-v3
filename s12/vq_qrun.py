"""s12 QUANTUM-ROLE -- the quantum arms.

stage2  EXACT-GRADIENT CVaR free-energy VQE (`core.quantum.run_cvar_vqe`): the alpha/T
        sweep, state entropy (so a collapse to the argmin is visible), point vs
        distribution readouts.  Expensive -- run on a stratified subset.
stage3  DEVICE-REALISTIC CVaR-VQE (`core.quantum.run_global_cvar_vqe`, lightning.qubit,
        SPSA, sampled CVaR) against classical searches under the project's OWN shared
        evaluation budget (`budget.BudgetedEnergyModel`: one charge per unique bitstring).
        This is the matched-budget comparison.

Usage:
  python -m s12.vq_qrun stage2 <family> <n_targets> [iters]
  python -m s12.vq_qrun stage3 <family> <n_targets> <budget1,budget2,...>
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
from s12 import vq_quantum as Qm
from s12.vq_build import FAMILIES
from s12.vq_run import native

ONEHOT = {"A2_8", "A3_5"}


def energy_of(fam, inst):
    """(E vector, n_qubits, decode-to-config map or None) for a family's encoding."""
    if fam in ONEHOT:
        H = V.onehot_ising(inst)
        E, B = V.energy_vector_onehot(inst, H)
        return E, H["N"], None, H
    E, cfg, N = V.energy_vector_index(inst)
    return E, N, cfg, None


def state_to_cfg(fam, inst, state, cfgmap, H):
    """Basis-state index -> assembly configuration (None if infeasible under one-hot)."""
    if cfgmap is not None:
        return cfgmap[int(state)]
    k, m = inst["k"], inst["m"]
    bits = ((int(state) >> np.arange(k * m - 1, -1, -1)) & 1).astype(int)
    out = np.zeros(k, int)
    for s in range(k):
        blk = bits[s * m:(s + 1) * m]
        if blk.sum() != 1:
            return None
        out[s] = int(np.argmax(blk))
    return out


def stratified(tg, n, k, minlen):
    """n targets spanning folds, lengths and the FAIL18 / other-108 split."""
    ok = [t for t in tg if V.segments(t["n"], k, minlen) is not None]
    f18 = [t for t in ok if t["pdb"] in I.FAIL18]
    oth = [t for t in ok if t["pdb"] not in I.FAIL18]
    nf = max(1, round(n * len(f18) / len(ok)))
    pick = [f18[i] for i in np.linspace(0, len(f18) - 1, min(nf, len(f18))).astype(int)]
    pick += [oth[i] for i in np.linspace(0, len(oth) - 1, n - len(pick)).astype(int)]
    return pick


def gibbs_at_entropy(E, H_target, lo=1e-4, hi=1e3, iters=80):
    """Boltzmann p ~ exp(-E/T) with T chosen so that H(p) equals `H_target` (in bits).

    THE decisive classical control for the CVaR distribution claim.  A free energy
    `<E> - T H` is minimised over the WHOLE simplex by exactly this distribution; the
    CVaR free energy the VQE minimises differs only in which part of the energy
    distribution is scored.  Matching entropy removes "the quantum state is just broader"
    as an explanation, so what is left is whether the SHAPE of the VQE's distribution
    carries anything a soft-min does not.
    """
    E = np.asarray(E, float)
    E = E - E.min()
    def ent(T):
        w = np.exp(-E / max(T, 1e-12))
        w = w / w.sum()
        nz = w[w > 0]
        return float(-(nz * np.log2(nz)).sum()), w
    for _ in range(iters):
        mid = math.sqrt(lo * hi)
        h, _ = ent(mid)
        if h < H_target:
            lo = mid
        else:
            hi = mid
    h, w = ent(math.sqrt(lo * hi))
    return w, h


def topk_states(E, k):
    return np.argsort(np.asarray(E, float), kind="stable")[:int(k)]


# ------------------------------------------------------------------------- stage 2
def stage2(fam, ntarg=12, iters=60, layers=3, alphas=(0.05, 0.1, 0.25, 0.5, 1.0),
           Ts=(0.0, 0.02, 0.1, 0.3), out=None):
    cfg = FAMILIES[fam]
    tg = stratified(I.targets(), ntarg, cfg["k"], cfg["minlen"])
    rows = []
    t00 = time.time()
    for q, t in enumerate(tg):
        if I.free_gb() < 1.5:
            print(f"  free_gb {I.free_gb():.2f} < 1.5 -- waiting", flush=True)
            while I.free_gb() < 1.5:
                time.sleep(20)
        inst = V.cached_instance(t, **cfg)
        E, N, cfgmap, H = energy_of(fam, inst)
        nat = native(t["pdb"])
        exh = C.a_exhaustive(inst["h"], inst["J"], inst["npairs"])
        Eall, cfgs = exh.pop("all_E"), exh.pop("all_cfg")
        Xall = V.assemble_batch(inst, cfgs)
        rr = I.kabsch_rmsd_batch(Xall, nat)
        opt = float(Eall.min())
        row = dict(pdb=t["pdb"], n=t["n"], fold=t["fold"], n_qubits=N, dim=int(len(E)),
                   opt=opt, E_min_register=float(E.min()),
                   oracle_best_rmsd=float(rr.min()),
                   exact_opt_rmsd=float(rr[int(np.argmin(Eall))]),
                   anchor_rmsd=float(I.ca_rmsd(inst["anchor"], nat)),
                   grad_audit=Qm.vqe_gradient_audit(E, N, layers=layers, alpha=0.2),
                   arms=[])
        for a in alphas:
            for T in Ts:
                r = Qm.run(E, alpha=a, T=T, n=N, layers=layers, iters=iters, seed=0)
                ro = Qm.readouts(E, r["p"])
                # point readout
                c_arg = state_to_cfg(fam, inst, ro["argmax_state"], cfgmap, H)
                c_sup = state_to_cfg(fam, inst, ro["best_in_support_state"], cfgmap, H)
                rm_arg = float(I.ca_rmsd(V.assemble(inst, c_arg), nat)) if c_arg is not None else None
                rm_sup = float(I.ca_rmsd(V.assemble(inst, c_sup), nat)) if c_sup is not None else None
                # DISTRIBUTION readout: p-weighted coordinate average over feasible states
                if cfgmap is not None:
                    Xw = np.tensordot(r["p"], Xall, axes=(0, 0))
                    feas_mass = 1.0
                else:
                    w = np.zeros(len(cfgs))
                    for st_, pv in enumerate(r["p"]):
                        if pv <= 0:
                            continue
                        cc = state_to_cfg(fam, inst, st_, cfgmap, H)
                        if cc is not None:
                            w[int(np.ravel_multi_index(tuple(cc), [inst["m"]] * inst["k"]))] += pv
                    feas_mass = float(w.sum())
                    Xw = np.tensordot(w / max(feas_mass, 1e-300), Xall, axes=(0, 0))
                rm_pw = float(I.ca_rmsd(Xw, nat))
                # ---- classical controls at MATCHED entropy -------------------------
                pg, hg = gibbs_at_entropy(E, r["entropy_bits"])
                if cfgmap is not None:
                    Xg = np.tensordot(pg, Xall, axes=(0, 0))
                else:
                    wg = np.zeros(len(cfgs))
                    for st_, pv in enumerate(pg):
                        if pv <= 1e-18:
                            continue
                        cc = state_to_cfg(fam, inst, st_, cfgmap, H)
                        if cc is not None:
                            wg[int(np.ravel_multi_index(tuple(cc), [inst["m"]] * inst["k"]))] += pv
                    Xg = np.tensordot(wg / max(wg.sum(), 1e-300), Xall, axes=(0, 0))
                rm_gibbs = float(I.ca_rmsd(Xg, nat))
                # coverage: best true RMSD among the 8 most probable / 8 lowest-energy states
                def _cover(pv_or_E, by_prob):
                    idx = (np.argsort(-pv_or_E)[:8] if by_prob else topk_states(pv_or_E, 8))
                    best = np.inf
                    for st_ in idx:
                        cc = state_to_cfg(fam, inst, int(st_), cfgmap, H)
                        if cc is None:
                            continue
                        j = int(np.ravel_multi_index(tuple(cc), [inst["m"]] * inst["k"]))
                        best = min(best, float(rr[j]))
                    return None if not np.isfinite(best) else best
                cov_vqe = _cover(r["p"], True)
                cov_gibbs = _cover(pg, True)
                cov_energy = _cover(E, False)
                row["arms"].append(dict(alpha=a, T=T, cvar=r["cvar"],
                                        rmsd_gibbs_matched=rm_gibbs,
                                        gibbs_entropy_bits=hg,
                                        cover8_vqe=cov_vqe, cover8_gibbs=cov_gibbs,
                                        cover8_lowest_energy=cov_energy,
                                        entropy_bits=r["entropy_bits"],
                                        max_entropy_bits=float(N),
                                        wall=r["wall"],
                                        argmax_energy=ro["argmax_energy"],
                                        best_in_support_energy=ro["best_in_support_energy"],
                                        support_size=ro["support_size"],
                                        p_top=ro["p_top"],
                                        gap_argmax=float(ro["argmax_energy"] - float(E.min())),
                                        gap_best_in_support=float(ro["best_in_support_energy"] - float(E.min())),
                                        feasible_mass=feas_mass,
                                        rmsd_argmax=rm_arg, rmsd_best_in_support=rm_sup,
                                        rmsd_p_weighted=rm_pw))
        rows.append(row)
        print(f"  [{fam} stage2] {q+1}/{len(tg)} {t['pdb']} {time.time()-t00:.0f}s "
              f"free={I.free_gb():.2f}", flush=True)
        I.write(out or f"vq_stage2_{fam}", dict(family=fam, cfg=cfg, iters=iters,
                                                layers=layers, rows=rows))
    print("wrote", out or f"vq_stage2_{fam}", flush=True)
    return rows


# ------------------------------------------------------------------------- stage 3
def stage3(fam, ntarg=None, budgets=(256, 512, 1024, 2048, 3200), layers=3, alpha=0.25,
           shots=16, restarts=1, seed=0, out=None):
    import warnings
    from core import quantum as Q
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    cfg = FAMILIES[fam]
    tg = [t for t in I.targets() if V.segments(t["n"], cfg["k"], cfg["minlen"]) is not None]
    if ntarg:
        tg = stratified(I.targets(), ntarg, cfg["k"], cfg["minlen"])
    rows = []
    t00 = time.time()
    for q, t in enumerate(tg):
        inst = V.cached_instance(t, **cfg)
        E, N, cfgmap, H = energy_of(fam, inst)
        nat = native(t["pdb"])
        exh = C.a_exhaustive(inst["h"], inst["J"], inst["npairs"])
        Eall, cfgs = exh.pop("all_E"), exh.pop("all_cfg")
        Xall = V.assemble_batch(inst, cfgs)
        rr = I.kabsch_rmsd_batch(Xall, nat)
        opt = float(E.min())
        row = dict(pdb=t["pdb"], n=t["n"], fold=t["fold"], n_qubits=N, dim=int(len(E)),
                   opt=opt, exact_opt_rmsd=float(rr[int(np.argmin(Eall))]),
                   oracle_best_rmsd=float(rr.min()),
                   anchor_rmsd=float(I.ca_rmsd(inst["anchor"], nat)), budgets={})
        for B in budgets:
            entry = {}
            ham = Qm.qubo_hamiltonian(E, N, eval_budget=B)
            t0 = time.perf_counter()
            try:
                res = Q.run_global_cvar_vqe(ham, layers=layers, alpha=alpha, shots=shots,
                                            restarts=restarts, seed=seed,
                                            optimizer="SPSA", device="lightning.qubit",
                                            final_shots=max(256, shots * 4), verbose=False)
                entry["vqe"] = dict(best_energy=float(res["best_seen_energy"]),
                                    readout_energy=float(res["vqe_energy"]),
                                    modal_energy=float(res["vqe_modal_energy"]),
                                    entropy_bits=float(res["distribution_entropy_bits"]),
                                    top1=float(res["distribution_top1_prob"]),
                                    spsa_iters=int(res["n_spsa_iterations_total"]),
                                    n_energy_evaluations=int(res["n_energy_evaluations"]),
                                    terminated_by=res["terminated_by"],
                                    wall=float(time.perf_counter() - t0))
            except Exception as e:                       # noqa: BLE001
                entry["vqe"] = dict(error=str(e), wall=float(time.perf_counter() - t0))
            for nm, fn in (("random", Qm.budgeted_random),
                           ("anneal", Qm.budgeted_anneal),
                           ("greedy_ls", Qm.budgeted_greedy_ls)):
                hm = Qm.qubo_hamiltonian(E, N, eval_budget=B)
                entry[nm] = fn(hm, seed=seed)
            row["budgets"][str(B)] = entry
        rows.append(row)
        print(f"  [{fam} stage3] {q+1}/{len(tg)} {t['pdb']} {time.time()-t00:.0f}s "
              f"free={I.free_gb():.2f}", flush=True)
        if q % 5 == 0 or q == len(tg) - 1:
            I.write(out or f"vq_stage3_{fam}", dict(family=fam, cfg=cfg, alpha=alpha,
                                                    shots=shots, layers=layers,
                                                    restarts=restarts, rows=rows))
    I.write(out or f"vq_stage3_{fam}", dict(family=fam, cfg=cfg, alpha=alpha, shots=shots,
                                            layers=layers, restarts=restarts, rows=rows))
    return rows


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "stage2":
        stage2(sys.argv[2], int(sys.argv[3]),
               iters=int(sys.argv[4]) if len(sys.argv) > 4 else 60)
    elif what == "stage3":
        nt = sys.argv[3]
        stage3(sys.argv[2], None if nt in ("all", "0") else int(nt),
               budgets=tuple(int(x) for x in sys.argv[4].split(",")))
