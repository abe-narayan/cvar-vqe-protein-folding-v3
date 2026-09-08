"""SPRINT 16 / QPHASE -- the objective-quality PHASE BOUNDARY for CVaR-VQE.

    python -m s16.qphase verify     # every pre-condition, before any claim
    python -m s16.qphase real       # T2: where the REAL objectives sit (rho, global + tail)
    python -m s16.qphase ladder     # T1: the boundary -- blend and noise, all controls
    python -m s16.qphase lr         # T3: is the CVaR baseline defect just over-stepping?
    python -m s16.qphase alpha      # T4: is alpha a tail parameter or a temperature?
    python -m s16.qphase doing      # T5: optimising / sampling / exploring / representing
    python -m s16.qphase report     # every table, from the checkpoints

Every mode checkpoints after each target, so an interrupted run is readable rather than
lost.  Every blended or noised arm is an ORACLE DIAGNOSTIC -- it reads the native to build
the quality knob -- and is labelled as such in every row it appears in.  The BOUNDARY is the
deliverable; no blended arm is ever a predictive result.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

from s16 import qphase_lib as L
from s14 import vqe_lib as V
from s14 import vqe_run as R
from s15 import seed as SD
from s12 import instrument as I

TARGETS = L.TARGETS19
TARGETS9 = L.TARGETS9
SEEDS = (0, 1, 2)
#: the mixing-weight rungs.  The REALISED rho is measured at every one of them and is the
#: only quality axis any result is reported on -- the docstring of `blend_objective` warns
#: that the map from mixing weight to rho is strongly non-linear, and it is: on this
#: instrument signal 0.2 buys +0.13 of rho on one target and +0.24 on another.
SIGNALS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0)
#: the independent second parameterisation.  sigma is on the uniformised truth, so the
#: realised rho depends on sigma alone and is nearly target-independent -- which is exactly
#: what makes it a clean cross-check of a boundary measured with the blend.
SIGMAS = (4.0, 2.0, 1.0, 0.5, 0.3, 0.2, 0.12, 0.06, 0.0)


def _t(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ============================================================================ VERIFY
def verify():
    """Every pre-condition this workstream's claims rest on, checked before any of them."""
    out = {}
    ins = L.inst("1CS9")

    # 1. the hoisted blend is the library call, bit for bit
    ub = V.uniformise(ins.hamil())
    d = []
    for s in (0.0, 0.37, 1.0):
        a = ins.blend(ub, s)
        b = V.blend_objective(ins.hamil(), ins.rmsd, s)
        d.append(float(np.abs(a - b).max()))
    out["blend_hoist_max_abs_diff"] = d
    assert max(d) == 0.0, d

    # 2. the hoisted noise is the library call
    r1 = np.random.default_rng(7); r2 = np.random.default_rng(7)
    a = ins.noisy(0.3, r1); b = V.noisy_truth(ins.rmsd, 0.3, r2)
    out["noise_hoist_max_abs_diff"] = float(np.abs(a - b).max())
    assert out["noise_hoist_max_abs_diff"] == 0.0

    # 3. signal=0 preserves the base ranking exactly, signal=1 the truth's
    out["rho_blend0_vs_base"] = float(V.spearman(ins.blend(ub, 0.0), ins.hamil()))
    out["rho_blend1_vs_truth"] = float(V.spearman(ins.blend(ub, 1.0), ins.rmsd))

    # 4. the rebuilt CA trace is the one the stored RMSD scored
    rng = np.random.default_rng(0)
    idx = rng.integers(0, ins.N, 256)
    got = I.kabsch_rmsd_batch(ins.ca(idx), ins.nat)
    out["ca_rebuild_max_abs_err"] = float(np.abs(np.asarray(got) - ins.rmsd[idx]).max())
    assert out["ca_rebuild_max_abs_err"] < 1e-6

    # 5. THE STRUCTURAL IDENTITY: a reweighting of a drawn set cannot beat its argmin
    E = ins.hamil()
    ud = L.untrained_draws(ins.n_qubits, 4096, 0)
    ti = L.tilt_samples(ud, E, 6.0, np.random.default_rng(0))
    out["tilt_samples_argmin_equals_untrained"] = bool(
        L.argmin_readout(E, ins.rmsd, ti)["best_e"] >=
        L.argmin_readout(E, ins.rmsd, ud)["best_e"] - 1e-12)
    out["tilt_samples_subset_of_untrained"] = bool(
        np.isin(np.unique(ti), np.unique(ud)).all())

    # 6. the exact untrained distribution matches the empirical draws
    p0 = L.untrained_probs(ins.n_qubits, 0)
    big = L.untrained_draws(ins.n_qubits, 8192, 0)
    c = np.bincount(big, minlength=ins.N) / big.size
    out["p0_vs_empirical_l1"] = float(np.abs(c - p0).sum())
    out["p0_sum"] = float(p0.sum())
    out["p0_entropy_bits"] = float(-(p0[p0 > 0] * np.log2(p0[p0 > 0])).sum())

    # 7. expected-distinct closed form against a simulation
    rng = np.random.default_rng(1)
    sim = np.mean([np.unique(rng.choice(ins.N, size=2048, p=p0)).size for _ in range(8)])
    out["expected_distinct_closed_form"] = L.expected_distinct(p0, 2048)
    out["expected_distinct_simulated"] = float(sim)

    # 8. the tie rule actually fires somewhere
    ties = int((ins.hamil() == ins.hamil().min()).sum())
    out["hamil_argmin_ties_1CS9"] = ties

    # 9. seeds are stable across processes (the `hash()` defect)
    out["stable_seed_check"] = int(SD.stable_seed("1CS9", 3, salt=L.SALT))

    # 10. Sprint 15's instrument constants still reproduce
    try:
        base = I.load_baseline() if hasattr(I, "load_baseline") else None
        out["instrument"] = "not re-derived here (cache read); see s15 verify"
    except Exception as ex:
        out["instrument_error"] = str(ex)

    L.ck("verify", "checks", out)
    print(json.dumps(out, indent=1))
    return out


# =========================================================== T2: THE REAL OBJECTIVES
def real():
    """Where every REAL objective sits on the quality axis, on the boundary's instrument.

    Global rho AND in-tail rho are both reported, because the tail is where the optimiser
    concentrates and Sprint 15 found `tail - bulk` negative for every objective including
    the ORACLE.  If the two axes disagree about which side of the boundary an objective is
    on, that disagreement is the result.

    AMBER obeys the standing binding rule (`amber_kind == 0 AND amber_idx != snap_index`)
    and prints its surviving per-target n.
    """
    done = L.ck_load("real")
    for pdb in TARGETS:
        if pdb in done:
            continue
        L.wait_for_cpu(tag=f"real {pdb}")
        ins = L.inst(pdb)
        row = {"n": ins.n, "N": int(ins.N), "fold": ins.fold,
               "rmsd_best_in_space": float(ins.rmsd.min()),
               "rmsd_mean_in_space": float(ins.rmsd.mean())}
        objs = {"hamil": ins.hamil(), "legacy": ins.legacy, "prior": ins.prior}
        for nm, E in objs.items():
            row[nm] = L.rho_profile(E, ins.rmsd, seed=L.seed_of(pdb + nm))
        m = (ins.amber_kind == 0) & (ins.amber_idx != ins.snap_index)
        idx = ins.amber_idx[m]
        if idx.size >= 64:
            row["amber_total"] = L.rho_profile(ins.amber_total[m], ins.rmsd[idx],
                                               seed=L.seed_of(pdb + "amber"))
            row["amber_total"]["stratum"] = "amber_kind==0 AND amber_idx!=snap_index"
            row["amber_total"]["n_rows"] = int(idx.size)
            row["amber_total"]["n_rows_before_mask"] = int(ins.amber_idx.size)
        # ORACLE reference point: the truth itself
        row["ORACLE_truth"] = L.rho_profile(ins.rmsd, ins.rmsd, seed=0)
        L.ck("real", pdb, row)
        _t(f"real {pdb}: hamil rho {row['hamil']['rho_global']:+.3f} "
           f"tail {row['hamil']['rho_tail_0.01']:+.3f} "
           f"argmin_pct {row['hamil']['argmin_pct']:.3f}")
    return L.ck_load("real")


def real_from_qrestraint():
    """The generative / selective objectives Sprint 15 already measured on THIS instrument.

    `s15/results/qrestraint.json` carries per-target `ordering_skill` for `E_ls_pred`,
    `E_ml_pred`, `E_ls_pool`, `E_combined`, `E_ORACLE_true`, `S14_disto_bayes`, `legacy`,
    `prior` and `amber_total` on the same nineteen targets, the same full registers and the
    same tail definition.  Re-deriving them would change nothing and cost hours, so they are
    read and folded into the placement table with their provenance attached.
    """
    p = os.path.join(L.ROOT, "s15", "results", "qrestraint.json")
    if not os.path.exists(p):
        return {}
    with open(p) as fh:
        d = json.load(fh)
    out = {}
    for pdb, t in d.get("targets", {}).items():
        if "_skill" in t:
            out[pdb] = t["_skill"]
    return out


def realarms(targets=TARGETS, seeds=SEEDS):
    """T2, the DIRECT form: run the full arm bundle on the REAL objectives themselves.

    Placing a real objective on the boundary by quoting its rho assumes the boundary is a
    function of rho.  This mode does not assume it: every real objective gets the identical
    arm table the ladder's rungs get, so its VQE-minus-control difference is MEASURED and
    can be checked against what its rho predicts.  Five objectives, all NATIVE-FREE:

      hamil            the deployed structural Hamiltonian (retrieval prior + distogram)
      E_ls_pred        the GENERATIVE distance-restraint residual
      S14_disto_bayes  the SAME distogram consumed SELECTIVELY (Bayes risk)
      legacy           the genuine Legacy/MJ energy
      prior            the retrieval torsion prior alone

    Each is rank-uniformised first, exactly as every ladder rung is, so the comparison is
    between objectives and not between marginals.
    """
    done = L.ck_load("realarms")
    for pdb in targets:
        if pdb in done:
            continue
        L.wait_for_cpu(tag=f"realarms {pdb}")
        V.wait_for_memory(1.2, "realarms")
        ins = L.inst(pdb)
        objs = {"hamil": ins.hamil(), "legacy": ins.legacy, "prior": ins.prior}
        objs.update(ins.restraint_objectives())
        rows = {}
        for nm, raw in objs.items():
            E = V.uniformise(raw)
            e_pct = L.ordinal_pct(E)
            prof = L.rho_profile(E, ins.rmsd, seed=L.seed_of(f"ra{pdb}{nm}"))
            cells = {str(s): L.run_arms(ins, E, s, e_pct=e_pct) for s in seeds}
            rows[nm] = {"profile": prof, "seeds": cells}
            _t(f"  {pdb} {nm}: rho {prof['rho_global']:+.3f} "
               f"tail {prof['rho_tail_0.01']:+.3f} argmin_pct {prof['argmin_pct']:.3f}")
        L.ck("realarms", pdb, rows)
    return L.ck_load("realarms")


# =============================================================== T1: THE PHASE BOUNDARY
def _ladder_cell(ins, E, seed, tag, e_pct):
    a = L.run_arms(ins, E, seed, e_pct=e_pct)
    return a


def ladder(family="blend", targets=TARGETS, seeds=SEEDS):
    """THE BOUNDARY.  Sweep objective quality, measure the realised rho at every rung, and
    run CVaR-VQE against every control on the identical objective array.

    ORACLE DIAGNOSTIC: every rung above signal=0 / sigma=inf reads the native to build the
    knob.  The boundary is the deliverable.
    """
    tag = f"ladder_{family}"
    done = L.ck_load(tag)
    rungs = SIGNALS if family == "blend" else SIGMAS
    for pdb in targets:
        if pdb in done:
            continue
        L.wait_for_cpu(tag=f"{tag} {pdb}")
        V.wait_for_memory(1.2, tag)
        ins = L.inst(pdb)
        ub = V.uniformise(ins.hamil()) if family == "blend" else None
        rows = {}
        t0 = time.time()
        for g in rungs:
            if family == "blend":
                E = ins.blend(ub, g)
            else:
                E = ins.noisy(g, SD.stable_rng(pdb, g, "noise", salt=L.SALT))
            prof = L.rho_profile(E, ins.rmsd, seed=L.seed_of(f"{pdb}{family}{g}"))
            e_pct = L.ordinal_pct(E)
            cells = {}
            for s in seeds:
                cells[str(s)] = L.run_arms(ins, E, s, e_pct=e_pct)
            rows[str(g)] = {"profile": prof, "seeds": cells}
            _t(f"  {pdb} {family}={g}: rho {prof['rho_global']:+.3f} "
               f"tail {prof['rho_tail_0.01']:+.3f} argmin_pct {prof['argmin_pct']:.3f}")
        L.ck(tag, pdb, rows)
        _t(f"{tag} {pdb} done in {time.time()-t0:.0f}s")
    return L.ck_load(tag)


def sched(targets=TARGETS, seeds=SEEDS, rungs=(0.0, 0.3, 0.6, 1.0)):
    """DOES THE VQE'S ADVANTAGE OVER A FIXED-TEMPERATURE TILT REDUCE TO ITS SCHEDULE?

    A VQE samples a SEQUENCE of distributions, broad early and narrow late; a single
    Boltzmann tilt does not.  `tilt_anneal` gives the classical reweighting exactly that
    schedule over exactly the VQE's own untrained support, and nothing else.  If it
    reproduces the VQE, then whatever the VQE has over a fixed tilt is annealing -- and the
    real annealer is already in the table doing it far better and far cheaper.
    """
    done = L.ck_load("sched")
    for pdb in targets:
        if pdb in done:
            continue
        L.wait_for_cpu(tag=f"sched {pdb}")
        ins = L.inst(pdb)
        ub = V.uniformise(ins.hamil())
        rows = {}
        for g in rungs:
            E = ins.blend(ub, g)
            e_pct = L.ordinal_pct(E)
            prof = L.rho_profile(E, ins.rmsd, seed=L.seed_of(f"sc{pdb}{g}"))
            cells = {str(s): L.run_arms(
                ins, E, s, e_pct=e_pct,
                arms=("vqe", "untrained", "tilt_samples", "tilt_exact", "tilt_anneal",
                      "anneal")) for s in seeds}
            rows[str(g)] = {"profile": prof, "seeds": cells}
        L.ck("sched", pdb, rows)
        _t(f"sched {pdb} done")
    return L.ck_load("sched")


# ======================================================= T3: THE WEAK-OPTIMISER TEST
LRS = (0.15, 0.075, 0.0375, 0.019, 0.0094)


def lr(targets=TARGETS, seeds=SEEDS):
    """Is the CVaR gradient-baseline defect simply OPTIMISER STRENGTH?

    Sprint 15 recorded that the biased `baseline="tail"` estimator gives -0.266 where the
    corrected `baseline="const"` gives -0.145, and that its gradient norm is 2.4x smaller,
    i.e. that it acts as a de-facto step-size reduction.  If that is the whole mechanism,
    then scaling the learning rate DOWN on the corrected estimator must reproduce it.

    Arms: `const` at five learning rates spanning 16x, and `tail` at the reference rate.
    Measured against distinct-configuration count, draw entropy, gradient norm and the two
    structural readouts, on the REAL objective and on two ORACLE-DIAGNOSTIC quality rungs so
    the answer is not read at one point of the quality axis.
    """
    done = L.ck_load("lr")
    for pdb in targets:
        if pdb in done:
            continue
        L.wait_for_cpu(tag=f"lr {pdb}")
        ins = L.inst(pdb)
        ub = V.uniformise(ins.hamil())
        rows = {}
        for g in (0.0, 1.0):
            E = ins.blend(ub, g)
            e_pct = L.ordinal_pct(E)
            prof = L.rho_profile(E, ins.rmsd, seed=L.seed_of(f"lr{pdb}{g}"))
            cells = {}
            for s in seeds:
                base = L.run_arms(ins, E, s, arms=("untrained",), e_pct=e_pct)
                cells[str(s)] = {"untrained": base["untrained"]}
                for x in LRS:
                    a = L.run_arms(ins, E, s, arms=("vqe",), lr=x, baseline="const",
                                   e_pct=e_pct)
                    cells[str(s)][f"const_lr{x}"] = a["vqe"]
                a = L.run_arms(ins, E, s, arms=("vqe",), lr=0.15, baseline="tail",
                               e_pct=e_pct)
                cells[str(s)]["tail_lr0.15"] = a["vqe"]
            rows[str(g)] = {"profile": prof, "seeds": cells}
        L.ck("lr", pdb, rows)
        _t(f"lr {pdb} done")
    return L.ck_load("lr")


# ================================================================== T4: WHAT ALPHA IS
ALPHAS = (1.0, 0.5, 0.25, 0.1, 0.05, 0.02)


def alpha(targets=TARGETS, seeds=SEEDS):
    """Is `alpha` a genuine CVaR tail parameter, or a hidden temperature/diversity knob?

    THE DECISIVE TEST.  For every alpha, take the VQE's realised draw entropy and build the
    ENTROPY-MATCHED classical Boltzmann tilt of the untrained circuit's OWN samples.  If a
    one-parameter classical temperature reproduces alpha's entire effect on the readout,
    then alpha is a temperature; if the VQE's alpha curve separates from the entropy-matched
    tilt's curve anywhere, alpha carries something the tail parameter alone supplies.

    Both readouts are carried through, because Sprint 15's whole story is that they come
    apart: the ARGMIN readout the pipeline uses, and the ENSEMBLE coordinate average.
    """
    done = L.ck_load("alpha")
    for pdb in targets:
        if pdb in done:
            continue
        L.wait_for_cpu(tag=f"alpha {pdb}")
        ins = L.inst(pdb)
        ub = V.uniformise(ins.hamil())
        rows = {}
        for g in (0.0, 1.0):
            E = ins.blend(ub, g)
            e_pct = L.ordinal_pct(E)
            cells = {}
            for s in seeds:
                c = L.run_arms(ins, E, s, arms=("untrained",), e_pct=e_pct)
                for a in ALPHAS:
                    r = L.run_arms(ins, E, s, alpha=a,
                                   arms=("vqe", "tilt_samples", "tilt_exact"), e_pct=e_pct)
                    c[f"vqe_a{a}"] = r["vqe"]
                    c[f"tiltmatch_a{a}"] = r["tilt_samples"]
                    c[f"tiltx_a{a}"] = r["tilt_exact"]
                cells[str(s)] = c
            rows[str(g)] = {"seeds": cells,
                            "profile": L.rho_profile(E, ins.rmsd,
                                                     seed=L.seed_of(f"al{pdb}{g}"))}
        L.ck("alpha", pdb, rows)
        _t(f"alpha {pdb} done")
    return L.ck_load("alpha")


# =========================================================== T5: WHAT THE VQE IS DOING
def doing(targets=TARGETS9, seeds=SEEDS):
    """Optimising, sampling, exploring, or representing?  Four axes, never substituted.

    OPTIMISING   -- does it reach the objective's certified optimum? (`objective_gap`,
                    `p(argmin)` under the exact final distribution)
    SAMPLING     -- does the distribution it learns put mass where the objective is low?
                    (`E_p[E]` percentile against the untrained circuit's)
    EXPLORING    -- how many distinct configurations does it visit per unit budget?
    REPRESENTING -- does the learned distribution enrich near-native mass, and does the
                    MODE (what the variational state actually learned) beat the untrained
                    mode?  This is the only readout for which "the VQE did something" is a
                    meaningful claim.

    `exact_dist=True`, so the whole 2^18 distribution is enumerated rather than sampled.
    Restricted to the nine n=9 targets: at n=10 an exact 2^20 `probs` per run is 8.6 s and
    buys nothing this question needs.
    """
    done = L.ck_load("doing")
    for pdb in targets:
        if pdb in done:
            continue
        L.wait_for_cpu(tag=f"doing {pdb}")
        ins = L.inst(pdb)
        ub = V.uniformise(ins.hamil())
        rows = {}
        for g in (0.0, 0.3, 0.6, 1.0):
            E = ins.blend(ub, g)
            e_pct = L.ordinal_pct(E)
            prof = L.rho_profile(E, ins.rmsd, seed=L.seed_of(f"do{pdb}{g}"))
            cells = {}
            for s in seeds:
                an, th0 = L.theta0(ins.n_qubits, s)
                p0 = L.untrained_probs(ins.n_qubits, s)
                r = R.run(E, ins.n_qubits, L.ALPHA, L.BUDGET, shots=L.SHOTS,
                          ansatz=L.ANSATZ, seed=s, rmsd=ins.rmsd,
                          bits_per_res=ins.bits_per_res, exact_dist=False,
                          baseline="const", keep_seen=True)
                p = np.asarray(an.probs(np.asarray(r["theta"], float)), float)
                p = np.maximum(p, 0.0); p = p / p.sum()
                tied = np.flatnonzero(E == E.min())
                q = e_pct
                hp = float(-(p[p > 0] * np.log2(p[p > 0])).sum())
                cells[str(s)] = {
                    "objective_gap": float(np.min(E[r["seen"]]) - E.min()),
                    "p_argmin": float(p[tied].sum()),
                    "p0_argmin": float(p0[tied].sum()),
                    "Ep_E_pct": float(p @ q), "Ep0_E_pct": float(p0 @ q),
                    "entropy_bits": hp,
                    "entropy0_bits": float(-(p0[p0 > 0] * np.log2(p0[p0 > 0])).sum()),
                    "eff_support": float(2.0 ** hp),
                    "n_distinct": int(np.unique(r["seen"]).size),
                    "mode_rmsd": float(ins.rmsd[int(np.argmax(p))]),
                    "mode0_rmsd": float(ins.rmsd[int(np.argmax(p0))]),
                    "mean_rmsd_under_p": float(p @ ins.rmsd),
                    "mean_rmsd_under_p0": float(p0 @ ins.rmsd),
                    "pmass_below_2.0": float(p[ins.rmsd < 2.0].sum()),
                    "p0mass_below_2.0": float(p0[ins.rmsd < 2.0].sum()),
                    "enrich_2.0_vs_untrained": float(
                        p[ins.rmsd < 2.0].sum() / max(p0[ins.rmsd < 2.0].sum(), 1e-12)),
                    "rmsd_returned": float(L.argmin_readout(E, ins.rmsd,
                                                            r["seen"])["rmsd_returned"]),
                    "rmsd_best_seen": float(ins.rmsd[r["seen"]].min()),
                    "rmsd_of_certified_optimum": float(ins.rmsd[tied].mean()),
                }
            rows[str(g)] = {"profile": prof, "seeds": cells}
        L.ck("doing", pdb, rows)
        _t(f"doing {pdb} done")
    return L.ck_load("doing")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "verify"
    _t(f"qphase {mode}; CPU {L.cpu_pct():.0f}%")
    if mode == "verify":
        verify()
    elif mode == "real":
        real()
    elif mode == "realarms":
        realarms()
    elif mode == "sched":
        sched()
    elif mode == "ladder":
        ladder("blend")
    elif mode == "ladder_noise":
        ladder("noise")
    elif mode == "lr":
        lr()
    elif mode == "alpha":
        alpha()
    elif mode == "doing":
        doing()
    elif mode == "report":
        from s16 import qphase_report as REP
        REP.main()
    else:
        raise SystemExit(f"unknown mode {mode!r}")
    _t(f"qphase {mode} complete")
