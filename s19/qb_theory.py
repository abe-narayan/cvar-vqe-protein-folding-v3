"""SPRINT 19 / AGENT B -- T1/T2: what CVaR can and cannot design, stated and VERIFIED.

This module exists to damage this lane's own hypothesis before the empirical arms are read.
It is pre-registered in `s19/PREREG_B.md` section 2 as prior (P-a).

T1  THE CVaR INDIFFERENCE THEOREM  (EXACT)
------------------------------------------
Let `E : Z -> R` with `E* = min E` attained on `Z*`.  For a distribution `p` on `Z` write
`CVaR_alpha(p)` for the mean of the lowest `alpha`-fraction of `E` under `p`.  Then

    (i)  CVaR_alpha(p) >= E*  for every p, with EQUALITY iff  p(Z*) >= alpha ;
    (ii) hence the minimiser is the FACE  { p : p(Z*) >= alpha },  on which p's behaviour on
         its remaining 1 - alpha mass is COMPLETELY UNCONSTRAINED.

**CVaR does not design a population.  It constrains an alpha-tail and is indifferent to
everything else.**  Any diversity, coverage or multimodality a CVaR-trained sampler exhibits
is a property of the variational family and of partial convergence -- not of the objective.

This generalises Sprint 17 section 9b (which proved the delta-collapse for the `mps2f`
family on the k=4 lattice, using the fact that the family contains every basis state) to
ANY family and ANY representation, continuous included -- and it does so without the
representability premise, which is exactly the premise that FAILS in this sprint's
continuous construction (a von Mises component with capped `kappa` cannot be a delta).

T1b  THE RESTRICTED-FAMILY FORM  (EXACT, and it is the one that binds here)
--------------------------------------------------------------------------
Inside a family `{p_theta}` the minimiser need not be degenerate, so (ii) does not apply
verbatim.  What still holds, and is what matters:

    CVaR_alpha(p_theta) is a functional of the LOWER alpha-TAIL OF THE PUSHFORWARD OF
    p_theta UNDER E, and of nothing else.

So the training signal is blind to every property of the candidate set that is not visible
in the lower tail of `E`: diversity, near-native coverage, the set mean, and RMSD itself.
Two parameters whose energy-pushforward tails agree are indistinguishable to CVaR however
different the structures they emit.

T2  THE NUMERICAL VERIFICATION
------------------------------
`t1` exhibits, on real targets, TWO distributions with IDENTICAL CVaR_alpha (to machine
precision) and maximally different diversity and realised RMSD -- the theorem made concrete.
`t2` samples the ACTUAL `mps2f` family at random parameters and measures the rank
correlation between CVaR and each downstream quantity, plus the spread of those quantities
INSIDE a CVaR level set.  The alignment axis is named (Spearman), as brief section 6 requires.

RUN:  python -m s19.qb_theory t1
      python -m s19.qb_theory t2
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402
from s17 import q_lib as QL               # noqa: E402
from core import quantum as Q             # noqa: E402
from s19 import qb_lib as L               # noqa: E402

ALPHAS = (0.05, 0.25, 0.50)


def _spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    den = np.linalg.norm(ra) * np.linalg.norm(rb)
    return float(ra @ rb / den) if den > 0 else float("nan")


def _cvar_weighted(e, p, alpha):
    """CVaR_alpha of the law `p` on the atoms with energies `e` -- exact, no sampling."""
    v, _q, _t = Q.cvar_from_probs(np.asarray(e, float), np.asarray(p, float), float(alpha))
    return float(v)


def t1(pdbs=None, M=4096, tag="t1"):
    """T1 -- two CVaR-OPTIMAL laws, identical objective, opposite populations."""
    from core import project as pj
    pdbs = pdbs or [t["pdb"] for t in I.targets()[:12]]
    L.gather_full()
    rows = []
    for pdb in pdbs:
        tgt = L.target(pdb)
        rng = SD.stable_rng(pdb, "t1", salt=L.SALT)
        b = (rng.random((M, tgt["n"])) < tgt["wmarg"][:, 1][None, :]).astype(np.int64)
        phi, psi = L.draw_from_basins(b, tgt["mu"], tgt["kap"], rng)
        o = L.new_obj(tgt, budget=M)
        e, CA = o.raw(phi, psi)
        rr = I.kabsch_rmsd_batch(CA, tgt["nat"])                      # ORACLE
        star = int(np.argmin(e))
        row = {"pdb": pdb, "n": tgt["n"], "fold": tgt["fold"], "M": M,
               "E_star": float(e[star]), "rmsd_of_E_argmin_ORACLE": float(rr[star]),
               "alphas": {}}
        for a in ALPHAS:
            # -- law A: the delta on the objective's argmin (minimum diversity)
            pA = np.zeros(M); pA[star] = 1.0
            # -- law B: alpha mass on the argmin, the rest spread UNIFORMLY (maximum diversity)
            pB = np.full(M, (1.0 - a) / (M - 1)); pB[star] = a
            # -- law C: alpha on the argmin, the rest on the WORST configs (adversarial)
            pC = np.zeros(M); pC[star] = a
            worst = np.argsort(e)[::-1][:max(1, M // 4)]
            pC[worst] = (1.0 - a) / len(worst)
            out = {}
            for nm, p in (("A_delta", pA), ("B_uniform_tail", pB), ("C_worst_tail", pC)):
                cv = _cvar_weighted(e, p, a)
                # the population this law actually emits, at m = 75 i.i.d. draws
                pick = rng.choice(M, size=75, p=p / p.sum())
                W = CA[pick]
                md = QL.md_plane(W, tgt["nat"])
                avg, _bm = I.coordinate_average(W)
                out[nm] = {"cvar": cv, "D": md["D"], "M_err_ORACLE": md["M"],
                           "set_best_ORACLE": float(rr[pick].min()),
                           "set_mean_ORACLE": float(rr[pick].mean()),
                           "avg_ORACLE": float(I.ca_rmsd(np.asarray(avg, float),
                                                         tgt["nat"]))}
            out["cvar_spread"] = float(max(v["cvar"] for v in out.values()
                                           if isinstance(v, dict))
                                       - min(v["cvar"] for v in out.values()
                                             if isinstance(v, dict)))
            row["alphas"][str(a)] = out
        rows.append(row)
        print(f"[t1] {pdb}  " + "  ".join(
            f"a={a}: cvar spread {row['alphas'][str(a)]['cvar_spread']:.3e}  "
            f"D {row['alphas'][str(a)]['A_delta']['D']:.2f}/"
            f"{row['alphas'][str(a)]['B_uniform_tail']['D']:.2f}  "
            f"avg {row['alphas'][str(a)]['A_delta']['avg_ORACLE']:.2f}/"
            f"{row['alphas'][str(a)]['B_uniform_tail']['avg_ORACLE']:.2f}"
            for a in ALPHAS), flush=True)
    with open(os.path.join(L.RESULTS, f"qb_{tag}.json"), "w") as fh:
        json.dump({"rows": rows, "complete": len(rows) == len(pdbs),
                   "n_expected": len(pdbs)}, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return rows


def t2(pdbs=None, K=256, shots=512, tag="t2"):
    """T1b -- the family-restricted form, measured on the ACTUAL mps2f family.

    `K` random parameter vectors; for each, the induced continuous law's sampled CVaR and
    the population it emits.  Reports Spearman(CVaR, .) and the spread of each downstream
    quantity INSIDE a narrow CVaR band, which is the operational content of "CVaR is a
    functional of the energy tail and of nothing else".
    """
    pdbs = pdbs or [t["pdb"] for t in I.targets()[:8]]
    L.gather_full()
    rows = []
    for pdb in pdbs:
        tgt = L.target(pdb)
        rng = SD.stable_rng(pdb, "t2", salt=L.SALT)
        an = Q.MPSAnsatz(tgt["n"], layers=2, final_ry=True, entangler="cnot")
        o = L.new_obj(tgt, budget=K * shots + 10)
        rec = {a: {"cvar": [], "D": [], "gen": [], "avg": [], "setmean": [], "ent": []}
               for a in ALPHAS}
        for _k in range(K):
            th = np.pi / 2 + rng.normal(0.0, 1.2, an.n_params())
            bits = an.sample(th, shots, rng)
            phi, psi = L.draw_from_basins(bits, tgt["mu"], tgt["kap"], rng)
            e, CA = o.raw(phi, psi)
            rr = I.kabsch_rmsd_batch(CA, tgt["nat"])                  # ORACLE
            allb = Q.all_bitstrings(tgt["n"]) if tgt["n"] <= 16 else None
            ent = float("nan")
            if allb is not None:
                p = np.exp(np.asarray(an.logp(th, allb), float))
                p = np.maximum(p, 0); p /= p.sum()
                pos = p[p > 0]
                ent = float(-(pos * np.log2(pos)).sum())
            for a in ALPHAS:
                cv, _q, _t = Q.cvar(e, a)
                sel = np.argsort(e, kind="stable")[:75]
                md = QL.md_plane(CA[sel], tgt["nat"])
                avg, _bm = I.coordinate_average(CA[sel])
                rec[a]["cvar"].append(float(cv))
                rec[a]["D"].append(md["D"])
                rec[a]["gen"].append(float(rr.min()))
                rec[a]["setmean"].append(float(rr[sel].mean()))
                rec[a]["avg"].append(float(I.ca_rmsd(np.asarray(avg, float), tgt["nat"])))
                rec[a]["ent"].append(ent)
        row = {"pdb": pdb, "n": tgt["n"], "fold": tgt["fold"], "K": K, "alphas": {}}
        for a in ALPHAS:
            d = {k: np.asarray(v, float) for k, v in rec[a].items()}
            ordr = np.argsort(d["cvar"])
            band = ordr[:max(8, K // 8)]              # the best 12.5% by CVaR
            row["alphas"][str(a)] = {
                "rho_cvar_D": _spearman(d["cvar"], d["D"]),
                "rho_cvar_gen_ORACLE": _spearman(d["cvar"], d["gen"]),
                "rho_cvar_avg_ORACLE": _spearman(d["cvar"], d["avg"]),
                "rho_cvar_setmean_ORACLE": _spearman(d["cvar"], d["setmean"]),
                "rho_cvar_entropy": _spearman(d["cvar"], d["ent"]),
                "band_cvar_range": [float(d["cvar"][band].min()),
                                    float(d["cvar"][band].max())],
                "band_rel_cvar_spread": float((d["cvar"][band].max()
                                               - d["cvar"][band].min())
                                              / max(1e-12, abs(d["cvar"].mean()))),
                "band_D_range": [float(d["D"][band].min()), float(d["D"][band].max())],
                "band_avg_range_ORACLE": [float(d["avg"][band].min()),
                                          float(d["avg"][band].max())],
                "band_gen_range_ORACLE": [float(d["gen"][band].min()),
                                          float(d["gen"][band].max())],
                "all_avg_range_ORACLE": [float(d["avg"].min()), float(d["avg"].max())],
            }
        rows.append(row)
        r025 = row["alphas"]["0.25"]
        print(f"[t2] {pdb}  rho(cvar,D)={r025['rho_cvar_D']:+.3f}  "
              f"rho(cvar,avg)={r025['rho_cvar_avg_ORACLE']:+.3f}  "
              f"band avg {r025['band_avg_range_ORACLE'][0]:.2f}-"
              f"{r025['band_avg_range_ORACLE'][1]:.2f} A at cvar spread "
              f"{100*r025['band_rel_cvar_spread']:.1f}%", flush=True)
    with open(os.path.join(L.RESULTS, f"qb_{tag}.json"), "w") as fh:
        json.dump({"rows": rows, "complete": len(rows) == len(pdbs),
                   "n_expected": len(pdbs)}, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return rows


def t3(pdbs=None, R=64, tag="t3"):
    """T3 -- HOW MUCH OF THE OBJECTIVE CAN THE LATENT EVEN SEE?

    Any quantum sampler over a continuous space acts through a DISCRETE latent.  Decompose

        Var(E) = Var_b( E[E | b] )   +   E_b( Var(E | b) )
                 -------- between --------  -------- within --------

    with `R` continuous draws per latent `b`.  The between term is the ENTIRE signal the
    circuit's distribution can move; the within term is noise from the sampler's own
    continuous law and is invisible to `p_theta`.  The same decomposition is run on the
    ORACLE CA-RMSD, because a latent that controls the objective but not the structure is
    a different failure from one that controls neither.

    This is a property of the CONSTRUCTION, not of quantumness -- the identical ceiling
    binds the CNOT-free and the classical-chain latents.  It is reported as such.
    """
    pdbs = pdbs or [t["pdb"] for t in I.targets()[:16]]
    L.gather_full()
    rows = []
    for pdb in pdbs:
        tgt = L.target(pdb)
        rng = SD.stable_rng(pdb, "t3", salt=L.SALT)
        n = tgt["n"]
        nb = min(1 << n, 256)
        B = (rng.random((nb, n)) < 0.5).astype(np.int64) if (1 << n) > 256 else \
            Q.all_bitstrings(n).astype(np.int64)
        o = L.new_obj(tgt, budget=nb * R + 10)
        me, ve, mr, vr = [], [], [], []
        for b in B:
            bb = np.repeat(b[None], R, 0)
            phi, psi = L.draw_from_basins(bb, tgt["mu"], tgt["kap"], rng)
            e, CA = o.raw(phi, psi)
            rr = I.kabsch_rmsd_batch(CA, tgt["nat"])                   # ORACLE
            me.append(e.mean()); ve.append(e.var())
            mr.append(rr.mean()); vr.append(rr.var())
        me = np.asarray(me); ve = np.asarray(ve)
        mr = np.asarray(mr); vr = np.asarray(vr)
        row = {"pdb": pdb, "n": n, "fold": tgt["fold"], "n_latents": int(nb), "R": R,
               "E_between_frac": float(me.var() / (me.var() + ve.mean())),
               "RMSD_between_frac_ORACLE": float(mr.var() / (mr.var() + vr.mean())),
               "E_best_latent_mean": float(me.min()),
               "E_worst_latent_mean": float(me.max()),
               "rmsd_of_best_E_latent_ORACLE": float(mr[int(np.argmin(me))]),
               "rmsd_best_latent_ORACLE": float(mr.min()),
               "rho_latentE_latentRMSD_ORACLE": _spearman(me, mr)}
        rows.append(row)
        print(f"[t3] {pdb}  between-latent share of Var(E) = "
              f"{row['E_between_frac']:.3f}   of Var(RMSD) = "
              f"{row['RMSD_between_frac_ORACLE']:.3f}   "
              f"rho(latent E, latent RMSD) = {row['rho_latentE_latentRMSD_ORACLE']:+.3f}",
              flush=True)
    with open(os.path.join(L.RESULTS, f"qb_{tag}.json"), "w") as fh:
        json.dump({"rows": rows, "complete": len(rows) == len(pdbs),
                   "n_expected": len(pdbs)}, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return rows


def t4(pdbs=None, tag="t4"):
    """T4 -- THE MODEL CLASS, and why no sampling separation is available at this depth.

    EXACT.  `core.quantum.MPSAnsatz(n, layers=L, entangler="cnot")` states in its own
    docstring, and its `chi` attribute confirms, that every interior bond is exactly
    `2**L` with NO truncation.  At the deployed `L = 2` the state is a bond-dimension-4
    MPS, so the measured law `p_theta(b) = |psi(b)|^2` is a matrix-product distribution of
    bond dimension at most 16 -- i.e. a hidden Markov model with at most 16 latent states,
    exactly normalisable and exactly samplable classically in O(n * chi^2).

    **There is therefore no sampling separation available from this ansatz at any budget,
    by construction and not by measurement.**  What the entangler can buy is a richer
    correlation structure over the latent than a product model has -- which is a modelling
    question, not a quantum one, and is exactly what the classical `c_chain` control tests.

    What this function MEASURES is whether the trained circuit uses even that capacity:
    the latent law's TOTAL CORRELATION (multi-information)

        C(b) = sum_i H(b_i) - H(b)      bits

    for the trained entangled circuit, its CNOT-free twin, and the untrained circuit.
    """
    from s19 import qb_arms as AR
    pdbs = pdbs or [t["pdb"] for t in I.targets()[:16]]
    L.gather_full()
    rows = []
    for pdb in pdbs:
        tgt = L.target(pdb)
        n = tgt["n"]
        if n > 16:
            continue
        allb = Q.all_bitstrings(n)
        row = {"pdb": pdb, "n": n, "fold": tgt["fold"], "arms": {}}
        for nm, kw in (("trained_cnot", dict(alpha=0.25, ansatz="mps2f", train=True)),
                       ("trained_nocnot", dict(alpha=0.25, ansatz="mps2fn", train=True)),
                       ("untrained_cnot", dict(alpha=0.25, ansatz="mps2f", train=False))):
            an = Q.MPSAnsatz(n, layers=2, final_ry=True,
                             entangler=("none" if kw["ansatz"].endswith("n") else "cnot"))
            # `AR.run_circuit` does not return theta, so the identical training loop is
            # replayed here on the identical seed stream.  `qb_verify_t4` checks that the
            # replay reproduces `AR.run_circuit`'s own diagnostics bit for bit.
            rng = SD.stable_rng(pdb, 0, kw["ansatz"], kw["alpha"], salt=L.SALT)
            th = np.pi / 2 + rng.normal(0.0, 0.8, an.n_params())
            opt = Q.Adam(an.n_params(), lr=L.LR)
            o2 = L.new_obj(tgt)
            for _t in range(o2.budget // L.SHOTS):
                bits = an.sample(th, L.SHOTS, rng)
                phi, psi = L.draw_from_basins(bits, tgt["mu"], tgt["kap"], rng)
                e = o2(phi, psi)
                if e.size < L.SHOTS:
                    break
                if kw["train"]:
                    g, _v = Q.cvar_gradient(an, th, bits, e, kw["alpha"], baseline="const")
                    th = opt.step(th, g)
            p = np.exp(np.asarray(an.logp(th, allb), float))
            p = np.maximum(p, 0); p /= p.sum()
            pos = p[p > 0]
            H = float(-(pos * np.log2(pos)).sum())
            Hm = 0.0
            for i in range(n):
                q1 = float(p[allb[:, i] == 1].sum())
                q1 = min(max(q1, 1e-15), 1 - 1e-15)
                Hm += -(q1 * np.log2(q1) + (1 - q1) * np.log2(1 - q1))
            ref = AR.run_circuit(tgt, L.new_obj(tgt), seed=0, **kw)
            row["arms"][nm] = {"chi": int(an.chi), "entropy_bits": H,
                               "sum_marginal_entropy_bits": Hm,
                               "total_correlation_bits": float(Hm - H),
                               "total_corr_per_residue": float((Hm - H) / n),
                               "replay_matches_run_circuit":
                                   bool(abs(H - ref.get("entropy_bits", np.nan)) < 1e-9)}
        rows.append(row)
        a = row["arms"]
        print(f"[t4] {pdb} n={n} chi={a['trained_cnot']['chi']}  total correlation (bits): "
              f"trained+CNOT {a['trained_cnot']['total_correlation_bits']:.4f}   "
              f"trained-CNOT {a['trained_nocnot']['total_correlation_bits']:.4f}   "
              f"untrained {a['untrained_cnot']['total_correlation_bits']:.4f}", flush=True)
    with open(os.path.join(L.RESULTS, f"qb_{tag}.json"), "w") as fh:
        json.dump({"rows": rows, "complete": True, "n_expected": len(rows)}, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return rows


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else "t1"
    {"t1": t1, "t2": t2, "t3": t3, "t4": t4}[m]()
