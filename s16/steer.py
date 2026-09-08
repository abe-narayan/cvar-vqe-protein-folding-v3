"""SPRINT 16, coordinator -- THE FLAGSHIP: can the two native-free ingredients be combined?

THE HYPOTHESIS, AND THE THEORETICAL WORRY THAT HAS TO BE TESTED BEFORE IT.

Sprint 15 left the largest unexplained ceiling in the programme: an ORACLE rotation of the fit's own
error out of the loud half, at FIXED MAGNITUDE, reaches 1.855 A (-1.822 [-2.037, -1.613], W/L 121/5),
and the best of 22 native-free arms captured about 4% of it. Two native-free ingredients were then
found separately and never combined:

    ingredient 1   `theta_fit - theta_pool`, the channel-disagreement direction
                   |cos| 0.390 with the true torsion error against a 0.157 null
    ingredient 2   the emitted structure's QUIET (low-response) Jacobian subspace
                   bottom-half overlap with the native's at cos^2 0.827 against a 0.499 null

THE WORRY, STATED BEFORE THE EXPERIMENT SO IT CANNOT BE RETROFITTED. Quiet torsion directions are
DEFINED as those that barely move the superposed CA trace. RMSD is a function of the superposed CA
trace. So motion along quiet directions should barely change RMSD -- by construction, not by
accident. And the restraint objective is a function of CA distances only, so quiet motion should
barely change the objective either.

If both hold, then the fit's error along quiet directions is simultaneously **unconstrained by the
data** and **irrelevant to the score**, and the ORACLE rotation is not removing error at all: it is
replacing a loud error with a quiet one of equal magnitude, which lands near the native because quiet
displacements do not move the trace. That operation needs `theta_native`. It is a COUNTERFACTUAL
about where the error points, not a reachable move from `theta_fit`.

If that is right, the useful question is not "can we rotate the error" but "**can we shrink its LOUD
component**" -- and the loud component is exactly the part the restraints can see, because loud means
it moves the trace, which moves the distances, which the distogram constrains. In that case the
alignment ceiling is the restraint-accuracy problem in different coordinates, and the flagship closes
with a mechanism rather than with more failed arms.

**That is a falsifiable claim and this module tests it rather than assuming it.**

WHAT IS MEASURED

  A. THE DECOMPOSITION.  Split the fit's true torsion error `e` into loud and quiet components
     against the emitted structure's own Jacobian, and report how much RMSD each carries. Also
     report how well ||e_loud|| is predicted by the NATIVE-FREE restraint residual. Slack there --
     loud error not explained by the residual -- is the only exploitable room, and its size is the
     number that decides whether the flagship has anything to work with.

  B. THE STEERING LADDER.  Step from `theta_fit` along a direction, with the step size chosen
     LEAVE-FOLD-OUT on the four training folds:

       fit                the plain restraint fit                            NATIVE-FREE, control
       step_full          along `theta_pool - theta_fit`, unprojected        NATIVE-FREE
       step_loud_r        the same, projected onto the top-r LOUD subspace   NATIVE-FREE  <- flagship
       step_quiet         the same, projected onto the QUIET complement      NATIVE-FREE, null
       ORACLE_step_true   along the true error                               ORACLE ceiling
       ORACLE_loud_true   the true error projected onto the loud subspace    ORACLE decomposition

  `step_quiet` is the control that decides whether the framing above is right. **If stepping along
  quiet directions moves RMSD appreciably, the definition of "quiet" is not doing what it claims and
  every alignment result in the record needs re-reading.** It should do nothing.

  `ORACLE_loud_true` versus `ORACLE_step_true` says how much of the reachable gain lives in the loud
  subspace at all -- which bounds what any projected native-free steer could ever recover.

Run:
    python -m s16.steer
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
    os.environ.setdefault(_v, "2")

from s12 import instrument as I              # noqa: E402
from s15 import align_lib as A               # noqa: E402
from s15 import distcal as C                 # noqa: E402
from s15 import distgeo as D                 # noqa: E402
from s15 import pooldist as P                # noqa: E402
from s15 import seed as SD                   # noqa: E402

RESULTS = os.path.join(ROOT, "s16", "results")
os.makedirs(RESULTS, exist_ok=True)

STEPS = (0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0)
RANKS = (2, 4, 6, 10, 13)


def native_torsions(pdb):
    """Native (phi, psi) from the cached universe. ORACLE -- evaluation and diagnostics only."""
    u = I.load_univ(pdb)
    from core import geometry as geo
    nat = np.asarray(u["nat_ca"], float)
    return nat


def pool_torsions(pdb, n, fold):
    """Retrieval-conditioned circular-mean torsions. NATIVE-FREE."""
    from s14 import retprior as R
    PHI, PSI, _ = R.windows(pdb, "top75")
    return R.circ_mean(PHI, axis=0), R.circ_mean(PSI, axis=0)


def subspaces(phi, psi):
    """Right singular vectors of the superposed-CA Jacobian, loud first.

    `||J e|| / sqrt(n)` is the RMSD to first order, so the singular values order torsion
    directions by exactly how much they move the score.
    """
    J, _Jr, CA = A.sup_jacobian(phi, psi)
    _U, s, Vt = np.linalg.svd(J, full_matrices=False)
    return J, Vt.T, s, CA          # V columns ordered by decreasing singular value


def rmsd_of(phi, psi, nat):
    return float(I.ca_rmsd(I.build_ca(phi, psi), nat))


def run(targets=None, n_start=4):
    from s14 import ladder as L
    tg = targets if targets is not None else I.targets()
    pdbs = [t["pdb"] for t in tg]
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    data = C.gather(tg)

    # leave-fold-out separation debias, exactly as every other module fits it
    deb = {}
    for f in sorted(set(folds)):
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(sp), sp))

    rows = []
    path = os.path.join(RESULTS, "steer.json")

    for c, t in enumerate(tg):
        p = t["pdb"]; d = data[p]
        i, j, n, fold = d["i"], d["j"], d["n"], d["fold"]
        sd = d["sd"]; nat = d["nat"]
        dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
        w = 1.0 / sd ** 2

        # ---- the fit, objective-selected among native-free starts
        starts = D.starts(p, d["seq"], n, fold, n_start, SD.stable_rng(p, "s16steer"))
        best = None
        ens_pairs = []
        for phi0, psi0, _tag in starts:
            ph, ps, f = A.fit(dhat, sd, i, j, phi0, psi0)
            ens_pairs.append((ph, ps))
            if best is None or f < best[0]:
                best = (f, ph, ps)
        _f, phi, psi = best
        th = np.concatenate([phi, psi])

        J, V, sv, CA = subspaces(phi, psi)
        nphi, npsi = _native_torsions(p, nat, d["seq"], fold)   # ORACLE, diagnostics only
        th_nat = np.concatenate([nphi, npsi])
        e = A.wrap(th_nat - th)                                    # ORACLE error vector

        ppl, psl = pool_torsions(p, n, fold)
        u = A.wrap(np.concatenate([ppl, psl]) - th)                # NATIVE-FREE direction

        # ---- A. the decomposition, per rank
        dec = {}
        for r in RANKS:
            Vl = V[:, :r]
            Pl = Vl @ Vl.T
            eL = Pl @ e; eQ = e - eL
            dec[r] = {
                "frac_norm_loud": float(np.linalg.norm(eL) / max(np.linalg.norm(e), 1e-12)),
                "frac_rmsd_loud": float(np.linalg.norm(J @ eL) / max(np.linalg.norm(J @ e), 1e-12)),
                "cos_u_e_loud": float(abs((Pl @ u) @ eL) /
                                      max(np.linalg.norm(Pl @ u) * np.linalg.norm(eL), 1e-12)),
            }

        #: CURVATURE. The loud/quiet split is a FIRST-ORDER statement at theta_fit, but the
        #: steering moves a finite distance and the Jacobian turns along the way. Compare the
        #: linear prediction ||J e|| / sqrt(n) against the RMSD actually realised at theta_fit,
        #: so no first-order share is ever read as a finite-step promise.
        lin = float(np.linalg.norm(J @ e) / np.sqrt(n))

        # the native-free restraint residual -- how visible is the loud error?
        dfit = np.sqrt(((CA[i] - CA[j]) ** 2).sum(1))
        resid = float(np.abs(dfit - dhat).mean())

        # ---- B. the steering ladder
        arms = {}
        m = 2 * n
        Pl4 = V[:, :4] @ V[:, :4].T
        #: STRICT quiet nulls. The complement of the top-4 still contains moderately loud
        #: directions, so it is not the null the framing needs; the bottom HALF and bottom
        #: QUARTER are. If either moves RMSD appreciably, "quiet" is not what the record says.
        Vq2 = V[:, m // 2:]; Pq2 = Vq2 @ Vq2.T
        Vq4 = V[:, 3 * m // 4:]; Pq4 = Vq4 @ Vq4.T
        #: a second NATIVE-FREE direction: toward the multi-start ensemble's own consensus,
        #: which is a different channel from the retrieval pool and is free here.
        #: circular mean of the multi-start ensemble's torsions -- a plain arithmetic mean is
        #: wrong on angles and would put the consensus in the wrong place near the branch cut.
        TH = np.array([np.concatenate([pp, qq]) for pp, qq in ens_pairs])
        ens = np.arctan2(np.sin(TH).mean(0), np.cos(TH).mean(0))
        uc = A.wrap(ens - th)
        dirs = {"full": u, "consensus": uc,
                "quiet_half": Pq2 @ u, "quiet_quarter": Pq4 @ u,
                "ORACLE_true": e, "ORACLE_loud4": Pl4 @ e,
                "ORACLE_quiet_half": Pq2 @ e}
        for r in RANKS:
            Vr = V[:, :r]
            dirs[f"loud{r}"] = (Vr @ Vr.T) @ u
            dirs[f"c_loud{r}"] = (Vr @ Vr.T) @ uc
        for name, vec in dirs.items():
            nv = float(np.linalg.norm(vec))
            if nv < 1e-12:
                arms[name] = {str(s): rmsd_of(phi, psi, nat) for s in STEPS}
                continue
            arms[name] = {}
            for s in STEPS:
                x = th + s * vec
                arms[name][str(s)] = rmsd_of(x[:n], x[n:], nat)

        rows.append({
            "pdb": p, "n": int(n), "fold": int(fold),
            "rmsd_fit": rmsd_of(phi, psi, nat),
            "resid": resid, "err_norm": float(np.linalg.norm(e)),
            "cos_u_e": float(abs(u @ e) / max(np.linalg.norm(u) * np.linalg.norm(e), 1e-12)),
            "sv_top4_share": float((sv[:4] ** 2).sum() / (sv ** 2).sum()),
            "lin_pred": lin,
            "decomp": dec, "arms": arms,
        })
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(path, "w") as fh:
                json.dump({"partial": rows, "n_done": c + 1}, fh)
            print(f"  {c+1}/{len(tg)} checkpointed", flush=True)

    # ---------------------------------------------------------------- aggregate
    inc = L.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    g = lambda k: np.asarray([r[k] for r in rows], float)          # noqa: E731
    base = g("rmsd_fit")

    out = {"n": len(rows), "incumbent": float(ref.mean()), "rows": rows,
           "steps": list(STEPS), "ranks": list(RANKS)}
    out["decomposition"] = {
        str(r): {k: float(np.mean([x["decomp"][r][k] for x in rows]))
                 for k in rows[0]["decomp"][RANKS[0]]} for r in RANKS}
    out["cos_u_e_mean"] = float(g("cos_u_e").mean())
    out["linear_vs_actual"] = {"lin_pred_mean": float(g("lin_pred").mean()),
                               "actual_mean": float(base.mean()),
                               "ratio": float(g("lin_pred").mean() / max(base.mean(), 1e-9))}
    out["sv_top4_share_mean"] = float(g("sv_top4_share").mean())

    # leave-fold-out step size per arm, chosen on the four training folds
    out["arms"] = {}
    for name in rows[0]["arms"]:
        picked, vals = {}, np.zeros(len(rows))
        for f in sorted(set(folds)):
            tr = folds != f
            bs = min(STEPS, key=lambda s: np.mean([rows[q]["arms"][name][str(s)]
                                                   for q in range(len(rows)) if tr[q]]))
            picked[int(f)] = bs
        for q, r in enumerate(rows):
            vals[q] = r["arms"][name][str(picked[int(folds[q])])]
        curve = {str(s): float(np.mean([r["arms"][name][str(s)] for r in rows])) for s in STEPS}
        out["arms"][name] = {
            "lfo_step": picked, "curve": curve,
            **I.summary(vals), "median": float(np.median(vals)),
            "vs_fit": I.paired(vals, base, folds=folds, names=pdbs),
            "vs_incumbent": I.paired(vals, ref, folds=folds, names=pdbs)}

    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s16_steer", out, n_expected=len(tg))

    print(f"\nn = {len(rows)}   incumbent {ref.mean():.3f}   plain fit {base.mean():.3f}")
    print(f"native-free |cos(u, e)| = {out['cos_u_e_mean']:.3f}   "
          f"top-4 singular share of the Jacobian = {out['sv_top4_share_mean']:.3f}\n")
    print("A. THE DECOMPOSITION -- how much of the error, and of the RMSD, is LOUD  (ORACLE)")
    print(f"  {'rank r':>8}{'||e_loud||/||e||':>20}{'RMSD share loud':>18}{'|cos(P u, e_loud)|':>22}")
    for r in RANKS:
        z = out["decomposition"][str(r)]
        print(f"  {r:>8}{z['frac_norm_loud']:>20.3f}{z['frac_rmsd_loud']:>18.3f}"
              f"{z['cos_u_e_loud']:>22.3f}")
    print("\nB. THE STEERING LADDER -- leave-fold-out step, target as the unit")
    print(f"  {'arm':<18}{'mean':>8}{'median':>8}{'vs plain fit':>26}{'step/fold':>22}")
    z = out["linear_vs_actual"]
    print(f"  first-order prediction ||J e||/sqrt(n) = {z['lin_pred_mean']:.3f} A against an "
          f"actual {z['actual_mean']:.3f} A  (ratio {z['ratio']:.3f})")
    print("  -- a ratio far from 1 means the error is too large for the linearisation.")
    for name in ("full", "consensus", "loud4", "loud10", "loud13", "c_loud10",
                 "quiet_half", "quiet_quarter",
                 "ORACLE_quiet_half", "ORACLE_loud4", "ORACLE_true"):
        if name not in out["arms"]:
            continue
        s = out["arms"][name]; v = s["vs_fit"]
        print(f"  {name:<18}{s['mean']:>8.3f}{s['median']:>8.3f}"
              f"  {v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]"
              f"{str(list(s['lfo_step'].values())):>22}")
    print("\n(`quiet` is the control that decides the framing: it should do NOTHING. If it moves"
          "\n RMSD appreciably, 'quiet' is not what the record says it is.)")
    return out


_NATCACHE = os.path.join(ROOT, "s16", "cache", "native_torsions.npz")


def _native_torsions(pdb, nat, seq, fold):
    """Native torsions consistent with the ideal-geometry builder. ORACLE, diagnostics only.

    The universes cache native CA coordinates, not native torsions, and the two are not in
    bijection -- an ideal-geometry chain cannot reproduce an arbitrary CA trace exactly. We
    therefore take the torsions of the CLOSEST ideal-geometry chain, which is what every other
    module in this project means by "the native torsions" and what the projection emits.
    Cached, because the projection costs seconds per target and this is a pure function of the
    native coordinates.
    """
    store = {}
    if os.path.exists(_NATCACHE):
        z = np.load(_NATCACHE)
        store = {k: z[k] for k in z.files}
    kphi, kpsi = pdb + "_phi", pdb + "_psi"
    if kphi in store:
        return store[kphi], store[kpsi]
    r = I.project(np.asarray(nat, float), seq, int(fold), lam=0.0, multi=True)
    store[kphi] = np.asarray(r["phi"], float)
    store[kpsi] = np.asarray(r["psi"], float)
    np.savez(_NATCACHE, **store)
    return store[kphi], store[kpsi]


if __name__ == "__main__":
    run()
