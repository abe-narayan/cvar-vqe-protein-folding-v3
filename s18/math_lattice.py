"""SPRINT 18 / MATH -- L1..L5: the lattice side of the bridge.

The brief's bridge is that first-order functional ANOVA under a product reference measure IS
the Walsh weight-<=1 projection when the measure is uniform on the enumerated lattice.  That
equivalence is the entire justification for the continuous version, so it is VERIFIED here,
not asserted.  Five things are measured, all EXACT (no sampling, no search: every target is
fully enumerated and every projection is closed form).

L1  EQUIVALENCE, QUBIT FACTORISATION.
    `anova1_discrete(E, n=n_q, k=2, uniform)` vs `walsh_truncate(E, n_q, 1)`.
    Claim: identical to machine precision.  If this fails the bridge is broken and F5 fires.

L2  EQUIVALENCE, RESIDUE FACTORISATION.
    `anova1_discrete(E, n=n, k=4, uniform)` vs `walsh_residue_additive(E, n, 2)`.
    Claim: the residue-additive ANOVA equals weight-0 + weight-1 + intra-residue weight-2.

L3  ORTHOGONALITY AND VARIANCE BUDGET.  `<E - E_le1, g> = 0` for every additive `g`, and
    `Var(E) = Var(E_le1) + Var(E_ge2)`.  Reported for both factorisations.  The two
    retained-variance fractions are the honest statement of how much of the "per-residue
    field" the strict truncation actually keeps.

L4  WHICH OBJECT WAS 2.411 A.  Certified argmin RMSD (ORACLE, post hoc) of:
      full `hamil_uniformised`, strict Walsh weight-<=1, residue-additive ANOVA,
      and -- because production never rank-transforms anything -- the same three on the RAW
      `hamil`.  Paired fold-clustered intervals, medians, W/L.

L5  ENCODING INVARIANCE, the decisive test of whether the strict object is physics.
    Relabelling the k = 4 torsion states of a residue is a pure bookkeeping change: it
    permutes the rows of a lookup table and changes no structure, no energy, and no ordering
    of configurations.  The residue-additive projection is EQUIVARIANT under it (proved and
    checked).  The strict Walsh weight-<=1 projection is NOT: it depends on which pairs of
    states share a leading bit.  So the strict object's argmin is sampled over random
    per-residue relabellings and its spread reported.  If 2.661 (the full objective) lies
    inside that spread, "degree-1 beats full" is a statement about a binary encoding.

RUN:  python -m s18.math_lattice run
      python -m s18.math_lattice mu       (L6: reference-measure sensitivity)
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

from s18 import math_lib as M

from s15 import seed as SD                       # noqa: E402
from s16 import qphase_lib as QP                 # noqa: E402
from s14 import vqe_lib as V                     # noqa: E402

TARGETS19 = QP.TARGETS19
N_PERM = 24            # random per-residue relabellings in L5
CONFIG = {"module": "s18/math_lattice.py", "targets": list(TARGETS19),
          "n_perm": N_PERM, "salt": M.SALT, "hamil_w": 0.25,
          "walsh_convention": "chi_S(x)=prod_{q in S}(-1)^{x_q}, qubit 0 = MSB"}


# --------------------------------------------------------------------------- L1-L4
def cell(pdb):
    ins = QP.inst(pdb)
    n, k, nq, N = ins.n, ins.k, ins.n_qubits, ins.N
    bpr = ins.bits_per_res
    E_raw = ins.hamil().astype(np.float64)
    E_uni = V.uniformise(E_raw)
    rmsd = ins.rmsd                                            # ORACLE, post hoc only
    out = {"pdb": pdb, "n": n, "k": k, "n_qubits": nq, "N": int(N),
           "fold": int(ins.fold), "bits_per_res": int(bpr),
           "rmsd_best_in_space_ORACLE": float(rmsd.min()),
           "rmsd_mean_in_space_ORACLE": float(rmsd.mean())}

    for tag, E in (("uniformised", E_uni), ("raw", E_raw)):
        # ---- L1: qubit-factorisation ANOVA == Walsh weight-<=1
        _, fq, A_q = M.anova1_discrete(E, nq, 2)               # uniform over {0,1}^nq
        W_q = M.walsh_truncate(E, nq, 1)
        d1 = float(np.max(np.abs(A_q - W_q)))
        scale = float(np.max(np.abs(E)) + 1e-300)

        # ---- L2: residue-factorisation ANOVA == weight<=1 + intra-residue weight-2
        E0, fr, A_r = M.anova1_discrete(E, n, k)               # uniform over {0..k-1}^n
        W_r = M.walsh_residue_additive(E, n, bpr)
        d2 = float(np.max(np.abs(A_r - W_r)))

        # ---- L3: orthogonality and the variance budget
        R_q = E - A_q
        R_r = E - A_r
        vE = float(E.var())
        # <residual, g> = 0 for every additive g: test on the basis {indicator(x_i = s)}
        S = M.states_of(N, n, k)
        orth_r = 0.0
        for i in range(n):
            for s in range(k):
                m = (S[:, i] == s)
                orth_r = max(orth_r, abs(float(R_r[m].mean() * m.mean())))
        bits = ((np.arange(N, dtype=np.int64)[:, None] >>
                 np.arange(nq - 1, -1, -1)[None, :]) & 1).astype(np.int8)
        orth_q = 0.0
        for q in range(nq):
            m = (bits[:, q] == 1)
            orth_q = max(orth_q, abs(float(R_q[m].mean() * m.mean())))

        row = {
            "L1_max_abs_diff_anova_vs_walsh_w1": d1,
            "L1_rel": d1 / scale,
            "L2_max_abs_diff_anova_vs_walsh_res": d2,
            "L2_rel": d2 / scale,
            "L3_var_frac_qubit_le1": float(A_q.var() / vE),
            "L3_var_frac_residue_le1": float(A_r.var() / vE),
            "L3_var_budget_resid_qubit": float(abs(A_q.var() + R_q.var() - vE) / vE),
            "L3_var_budget_resid_residue": float(abs(A_r.var() + R_r.var() - vE) / vE),
            "L3_max_orth_residual_qubit": orth_q,
            "L3_max_orth_residual_residue": orth_r,
            # ---- L4: where each object's certified argmin sits (ORACLE)
            "argmin_rmsd_full_ORACLE": float(rmsd[E == E.min()].mean()),
            "argmin_rmsd_walsh_le1_ORACLE": float(rmsd[A_q == A_q.min()].mean()),
            "argmin_rmsd_residue_add_ORACLE": float(rmsd[A_r == A_r.min()].mean()),
            "argmin_ties_full": int((E == E.min()).sum()),
            "argmin_ties_walsh_le1": int((A_q == A_q.min()).sum()),
            "argmin_ties_residue_add": int((A_r == A_r.min()).sum()),
            "rho_full_rmsd_ORACLE": M.spearman(E, rmsd),
            "rho_walsh_le1_rmsd_ORACLE": M.spearman(A_q, rmsd),
            "rho_residue_add_rmsd_ORACLE": M.spearman(A_r, rmsd),
            "rho_walsh_le1_with_full": M.spearman(A_q, E),
            "rho_residue_add_with_full": M.spearman(A_r, E),
            "coef_residue_field": fr.tolist(),      # persisted objective coefficients
            "E0": float(E0),
        }
        out[tag] = row

    # ---- L5: encoding invariance of the strict object (on the uniformised objective,
    #          which is the one Sprint 17 truncated)
    out["L5"] = perm_test(ins, E_uni, rmsd)
    return out


def perm_test(ins, E, rmsd, n_perm=N_PERM):
    """Relabel the k torsion states per residue at random; redo both projections.

    The permutation acts on the configuration index only.  `E_pi(x) = E(pi(x))` is the SAME
    function of the physical conformation -- the lattice point `pi(x)` is the conformation
    `x` used to be.  So any projection that is a property of the physics must give the same
    argmin CONFORMATION.  The residue-additive projection does; the strict Walsh one does
    not, because the qubit factorisation is not preserved by a permutation of 4 states.
    """
    n, k, nq, N = ins.n, ins.k, ins.n_qubits, ins.N
    bpr = ins.bits_per_res
    pw = k ** np.arange(n - 1, -1, -1)
    S = M.states_of(N, n, k)
    a_walsh, a_res, ident_walsh, ident_res = [], [], [], []
    for t in range(n_perm):
        rng = SD.stable_rng(ins.pdb, t, "perm", salt=M.SALT)
        P = np.stack([rng.permutation(k) for _ in range(n)])          # (n, k)
        Sp = P[np.arange(n)[None, :], S]                              # relabelled states
        jp = (Sp * pw[None, :]).sum(1)                                # new index of each x
        # E_pi(x) = E(pi(x)); rmsd_pi(x) = rmsd(pi(x)) -- same conformation, new label
        Ep = E[jp]
        rp = rmsd[jp]
        _, _, Aq = M.anova1_discrete(Ep, nq, 2)
        _, _, Ar = M.anova1_discrete(Ep, n, k)
        a_walsh.append(float(rp[Aq == Aq.min()].mean()))
        a_res.append(float(rp[Ar == Ar.min()].mean()))
        ident_walsh.append(int((Aq == Aq.min()).sum()))
        ident_res.append(int((Ar == Ar.min()).sum()))
    # the unpermuted values, for reference
    _, _, Aq0 = M.anova1_discrete(E, nq, 2)
    _, _, Ar0 = M.anova1_discrete(E, n, k)
    base_w = float(rmsd[Aq0 == Aq0.min()].mean())
    base_r = float(rmsd[Ar0 == Ar0.min()].mean())
    return {"n_perm": int(n_perm),
            "walsh_le1_argmin_rmsd_ORACLE": a_walsh,
            "residue_add_argmin_rmsd_ORACLE": a_res,
            "walsh_le1_mean": float(np.mean(a_walsh)),
            "walsh_le1_sd": float(np.std(a_walsh)),
            "walsh_le1_min": float(np.min(a_walsh)),
            "walsh_le1_max": float(np.max(a_walsh)),
            "residue_add_mean": float(np.mean(a_res)),
            "residue_add_sd": float(np.std(a_res)),
            "residue_add_invariant": bool(np.allclose(a_res, base_r, atol=1e-9)),
            "unpermuted_walsh_le1": base_w,
            "unpermuted_residue_add": base_r,
            "argmin_ties_walsh": ident_walsh, "argmin_ties_res": ident_res}


# ------------------------------------------------------------------------- L6: mu
def mu_candidates(ins):
    """Three NATIVE-FREE product reference measures on the k-state lattice.

    uniform     -- the lattice's own counting measure; the only one under which the qubit
                   factorisation is even defined (see below).
    pool        -- the retrieval pool's empirical per-residue torsion marginal,
                   `s14.retprior.state_prior(pdb, 'top75', k)['P']`.  Target-conditioned and
                   native-free: it is built from retrieved fragments only.
    rama        -- the shipped leave-fold-out Ramachandran log-density of the residue's
                   class, evaluated at the k lattice torsions and normalised.  Generic,
                   sequence-conditioned, native-free.

    IMPORTANT AND STRUCTURAL: `pool` and `rama` are product measures over RESIDUES but they
    do not factorise over QUBITS (a general distribution on 4 states is not a product of two
    Bernoullis).  So under either of them the strict Walsh weight-<=1 object is NOT DEFINED
    as an orthogonal projection at all.  That is not a technicality; it means the strict
    object cannot be carried to any informative reference measure.
    """
    from s14 import retprior as RP
    from core import project as PJ
    n, k = ins.n, ins.k
    out = {"uniform": np.full((n, k), 1.0 / k)}
    try:
        P = np.asarray(RP.state_prior(ins.pdb, arm="top75", k=k)["P"], float)
        P = np.clip(P, 1e-9, None)
        out["pool"] = P / P.sum(1, keepdims=True)
    except Exception as e:                                          # pragma: no cover
        out["pool_error"] = str(e)
    try:
        T = PJ.logp_tables("rama")[int(ins.fold)]
        cls = PJ.res_classes(ins.seq)
        L = T[cls]                                                  # (n, RB, RB)
        lp = PJ._bilinear(L, ins.PHI.T, ins.PSI.T)                  # (k, n)
        lp = np.asarray(lp, float).T                                # (n, k)
        w = np.exp(lp - lp.max(1, keepdims=True))
        out["rama"] = w / w.sum(1, keepdims=True)
    except Exception as e:                                          # pragma: no cover
        out["rama_error"] = str(e)
    return out


def mu_cell(pdb):
    ins = QP.inst(pdb)
    n, k, nq = ins.n, ins.k, ins.n_qubits
    E_uni = V.uniformise(ins.hamil().astype(np.float64))
    rmsd = ins.rmsd
    mus = mu_candidates(ins)
    out = {"pdb": pdb, "n": n, "fold": int(ins.fold)}
    for name in ("uniform", "pool", "rama"):
        if name not in mus:
            out[name] = {"error": mus.get(name + "_error", "missing")}
            continue
        p = mus[name]
        E0, f, A = M.anova1_discrete(E_uni, n, k, p)
        R = E_uni - A
        # variance fractions here are mu-WEIGHTED (the projection is mu-orthogonal)
        w = np.prod(p[np.arange(n)[None, :], M.states_of(ins.N, n, k)], axis=1)
        w = w / w.sum()
        mE = float(w @ E_uni)
        vE = float(w @ (E_uni - mE) ** 2)
        vA = float(w @ (A - float(w @ A)) ** 2)
        vR = float(w @ (R - float(w @ R)) ** 2)
        out[name] = {
            "var_frac_residue_le1_mu_weighted": vA / max(vE, 1e-300),
            "var_budget_resid": abs(vA + vR - vE) / max(vE, 1e-300),
            "argmin_rmsd_residue_add_ORACLE": float(rmsd[A == A.min()].mean()),
            "argmin_ties": int((A == A.min()).sum()),
            "rho_with_full": M.spearman(A, E_uni),
            "rho_with_rmsd_ORACLE": M.spearman(A, rmsd),
            "coef_residue_field": f.tolist(),
            "E0": float(E0),
            "mu": p.tolist(),
            "qubit_factorises": bool(_factorises(p, ins.bits_per_res)),
        }
    out["argmin_rmsd_full_ORACLE"] = float(rmsd[E_uni == E_uni.min()].mean())
    return out


def _factorises(p, bpr, tol=1e-9):
    """Does each residue's k-state law factorise into `bpr` independent bits?"""
    n, k = p.shape
    if bpr == 1:
        return True
    for i in range(n):
        q = p[i].reshape((2,) * bpr)
        m = q
        marg = []
        for a in range(bpr):
            ax = tuple(b for b in range(bpr) if b != a)
            marg.append(q.sum(axis=ax))
        prod = marg[0]
        for a in range(1, bpr):
            prod = np.multiply.outer(prod, marg[a])
        if not np.allclose(m, prod, atol=tol):
            return False
    return True


# ------------------------------------------------------------------------- driver
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    tg = sys.argv[2].split(",") if len(sys.argv) > 2 else list(TARGETS19)
    if mode == "run":
        done = M.ck_load("lattice")
        for pdb in tg:
            if f"cell_{pdb}" in done:
                print("  skip", pdb, flush=True)
                continue
            t0 = time.time()
            r = cell(pdb)
            M.ck("lattice", f"cell_{pdb}", r)
            done[f"cell_{pdb}"] = r
            u = r["uniformised"]
            print(f"  {pdb} {time.time()-t0:.0f}s  L1 {u['L1_rel']:.2e}  "
                  f"L2 {u['L2_rel']:.2e}  varfrac q {u['L3_var_frac_qubit_le1']:.3f} "
                  f"res {u['L3_var_frac_residue_le1']:.3f}  argmin "
                  f"full {u['argmin_rmsd_full_ORACLE']:.3f} "
                  f"w1 {u['argmin_rmsd_walsh_le1_ORACLE']:.3f} "
                  f"res {u['argmin_rmsd_residue_add_ORACLE']:.3f}", flush=True)
        if all(f"cell_{p}" in M.ck_load("lattice") for p in TARGETS19):
            print("  sealing:", M.seal("lattice", CONFIG))
    elif mode == "mu":
        done = M.ck_load("mu")
        for pdb in tg:
            if f"mu_{pdb}" in done:
                print("  skip", pdb, flush=True)
                continue
            t0 = time.time()
            r = mu_cell(pdb)
            M.ck("mu", f"mu_{pdb}", r)
            print(f"  {pdb} {time.time()-t0:.0f}s  " + "  ".join(
                f"{k}:{r[k].get('argmin_rmsd_residue_add_ORACLE', float('nan')):.3f}"
                f"(v{r[k].get('var_frac_residue_le1_mu_weighted', float('nan')):.2f})"
                for k in ("uniform", "pool", "rama")), flush=True)
        if all(f"mu_{p}" in M.ck_load("mu") for p in TARGETS19):
            print("  sealing:", M.seal("mu", CONFIG))
    else:
        raise SystemExit(mode)


if __name__ == "__main__":
    main()
