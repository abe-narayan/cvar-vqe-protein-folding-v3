"""SPRINT 15, ALIGN, TASKS 2 and 3 -- CAN ALIGNMENT BE ENGINEERED RATHER THAN INHERITED?

Every well-aligned emitter the INFO workstream measured is aligned because it is a REAL
FRAGMENT.  Nothing says a generator can be STEERED into the RMSD-quiet subspace.  This module
tests four concrete interventions against their stated nulls, and -- the part that decides
whether "alignment" is a lever or only a description -- measures BOTH axes for every arm:

    RAW ACCURACY   torsion RMS error (deg), emitted-distance MAE, restraint residual
    ALIGNMENT      ||J e|| / (||e|| s_rms), with J at the native AND at the emitted torsions
    OUTCOME        full-chain CA-RMSD, the frozen definition

If RMSD improves only where raw error also improves, alignment is a DESCRIPTION, not a lever,
and this file must say so.

THE ARMS.  Every arm is the identical fit -- identical starts (`s15.distgeo.starts`, stable
seed), identical exact analytic O(n) torsion gradient, identical OBJECTIVE-ONLY multi-start
selection, identical leave-fold-out separation-debiased restraints -- differing only in the
per-restraint weight, the loss, or the subspace it is allowed to move in.

  squared        w = 1/sd^2, squared loss                        THE CONTROL
  unit_w         w = 1, squared loss                             the brief's second null
  term<b>_t<t>   terminal restraints down-weighted by t, for pairs touching the outer b
                 residues.  Motivated by INFO C.2: at coverage 0.7 terminal gaps emit 1.716 A
                 and mid-chain gaps 3.505 A, a 1.79 A spread.  If terminal error is cheap,
                 do not spend fit accuracy buying it down.
  jacw_a<a>      w *= (||d d_p/d theta|| / median)^a.  Crude response weighting: a < 0
                 down-weights the restraints the structure responds to most strongly, which
                 pushes the residual error into directions that move the structure less.
  damp_k<k>      w /= (1 + k * (a_p / median a_p)^2), a_p = ||J M^{-1} g_p|| = the coordinate
                 damage per unit error in restraint p, with M = sum_p w_p g_p g_p'.  This is
                 the principled version of the same idea: it prices a restraint by how much
                 DAMAGE its error does, not by how strongly the structure responds to it.
  rank_f<f>      solve restricted to theta_start + span(top-r right singular vectors of J),
                 r = f * 2n.  The brief's "projection onto the well-conditioned subspace".
  cauchy_s<s>,
  welsch_s<s>,
  gnc_cauchy     the ROBUST arms from `s15/robust.py`, re-run here on the full 126 so that
                 robustness and alignment are measured on ONE start-matched instrument.
                 K7: median z = 1.336, max z = 7.237 -- the violations are concentrated.
  <combos>       cauchy x terminal, cauchy x damp, terminal x damp.  If robustness and
                 alignment are the SAME lever the combination is sub-additive; if they are
                 two levers it is additive.  That is the test, and it is pre-declared here.

  ORACLE_*       ceilings and diagnostics ONLY, never headline:
                 `ORACLE_kill_loud_<f>`   take the `squared` solution and delete the component
                                          of its TRUE torsion error lying in the loud
                                          (top-f) subspace.  This is the most alignment
                                          engineering could ever buy at fixed |error|.
                 `ORACLE_term_native`     set the outer residues' torsions to native.

HYPERPARAMETERS.  Every sweep value is run on every target, and the reported arm is assembled
LEAVE-FOLD-OUT afterwards: for held-out fold f the value is the one with the best mean RMSD on
the OTHER FOUR folds.  The in-fold best (`*_INFOLD`) is reported beside it as the tuning
ceiling, labelled ORACLE, exactly as `s15/distcal.py` does.  Nothing is tuned on the fold it
is reported on.

    python -m s15.align_fit            # ~90 min on this machine, checkpoints every 10 targets
"""
from __future__ import annotations

import json
import os
import sys
import time

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
from s15 import info_lib as L                # noqa: E402
from s15 import robust as RB                 # noqa: E402
from s15 import seed as SD                   # noqa: E402
import peptide_db as pdb                     # noqa: E402

RESULTS = os.path.join(ROOT, "s15", "results")
PATH = os.path.join(RESULTS, "align_fit.json")

#: sweep grids.  Trimmed from a wider first draft purely for wall-clock on a machine at 98%
#: CPU; the values kept span the sign of every effect, which is what the test needs.
TERM = ((2, 0.0), (2, 0.25), (1, 0.0))
JACW = (-0.5, 0.5)
DAMP = (0.5, 2.0)
RANKF = (0.5, 0.75)
CAUCHY = (1.0, 2.0, 3.0)
WELSCH = ()
KILLF = (0.25, 0.5)


def term_weight(n, i, j, b, tau):
    """Multiplicative weight: `tau` for any restraint touching the outer `b` residues."""
    edge = (i < b) | (j >= n - b)
    return np.where(edge, tau, 1.0)


def damage_weight(Jsup, Gp, w0, k):
    """`1 / (1 + k (a_p/median a_p)^2)` with `a_p = ||J M^{-1} g_p||`.  NATIVE-FREE."""
    M = (Gp * w0[:, None]).T @ Gp
    M = M + 1e-6 * np.trace(M) / M.shape[0] * np.eye(M.shape[0])
    X = np.linalg.solve(M, Gp.T)                       # (2n, P)
    a = np.linalg.norm(Jsup @ X, axis=0)               # (P,)
    med = float(np.median(a))
    if med <= 0:
        return np.ones(len(a))
    return 1.0 / (1.0 + k * (a / med) ** 2)


def run(targets=None, n_start=6, limit=None):
    tg = targets if targets is not None else I.targets()
    if limit:
        tg = tg[:limit]
    pdbs = [t["pdb"] for t in tg]
    data = C.gather(tg)

    # ---- leave-fold-out separation debias, exactly as s15/robust.py
    folds = sorted({data[p]["fold"] for p in pdbs})
    deb = {}
    for f in folds:
        train = [p for p in pdbs if data[p]["fold"] != f]
        fn, _ = C.fit_correction(data, train, "sep")
        deb[f] = (lambda sep, _fn=fn: -_fn(np.zeros_like(sep), sep))

    rows = []
    t_start = time.time()
    for c, t in enumerate(tg):
        p = t["pdb"]
        d = data[p]
        n = d["n"]
        i, j = d["i"], d["j"]
        sd = d["sd"]
        dhat = np.maximum(d["dhat"] - deb[d["fold"]](d["sep"]), 2.0)
        nat = d["nat"]
        dtrue = d["dtrue"]
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        w0 = 1.0 / (sd * sd)
        S = D.starts(p, d["seq"], n, d["fold"], n_start, SD.stable_rng(p))

        row = {"pdb": p, "n": n, "fold": d["fold"], "arms": {}}

        def emit(name, phi, psi):
            row["arms"][name] = A.score_axes(phi, psi, nphi, npsi, nat, dhat, sd, i, j, dtrue)

        def best_over_starts(fn):
            """`fn(phi0, psi0) -> (phi, psi, f)`; select on the OBJECTIVE only."""
            bb = None
            for phi0, psi0, _tag in S:
                ph, ps, fv = fn(phi0, psi0)
                if bb is None or fv < bb[0]:
                    bb = (fv, ph, ps)
            return bb[1], bb[2]

        # ---------------------------------------------------------------- controls
        sq_per_start = []
        bb = None
        for phi0, psi0, _tag in S:
            ph, ps, fv = A.fit(dhat, sd, i, j, phi0, psi0, wpair=w0)
            sq_per_start.append((ph, ps, fv))
            if bb is None or fv < bb[2]:
                bb = (ph, ps, fv)
        sq_phi, sq_psi = bb[0], bb[1]
        emit("squared", sq_phi, sq_psi)
        emit("unit_w", *best_over_starts(
            lambda a, b: A.fit(dhat, sd, i, j, a, b, wpair=np.ones(len(dhat)))))

        # ---------------------------------------------------------------- terminal relax
        for b, tau in TERM:
            wt = w0 * term_weight(n, i, j, b, tau)
            emit(f"term{b}_t{tau}", *best_over_starts(
                lambda a, bb_, _w=wt: A.fit(dhat, sd, i, j, a, bb_, wpair=_w)))

        # -------- stage-1 geometry, shared by every Jacobian-aware arm (NATIVE-FREE)
        Jsup, Jraw, CA1 = A.sup_jacobian(sq_phi, sq_psi)
        Gp, _d1 = A.pair_jacobian(Jraw, CA1, i, j)
        resp = np.linalg.norm(Gp, axis=1)
        rmed = float(np.median(resp))
        s_sv, V_sv = A.spectrum(Jsup)

        for al in JACW:
            wt = w0 * (np.maximum(resp, 1e-9) / max(rmed, 1e-9)) ** al
            emit(f"jacw_a{al}", *best_over_starts(
                lambda a, bb_, _w=wt: A.fit(dhat, sd, i, j, a, bb_, wpair=_w)))

        for k in DAMP:
            wt = w0 * damage_weight(Jsup, Gp, w0, k)
            emit(f"damp_k{k}", *best_over_starts(
                lambda a, bb_, _w=wt: A.fit(dhat, sd, i, j, a, bb_, wpair=_w)))

        for f in RANKF:
            r = max(1, int(round(f * 2 * n)))
            B = V_sv[:, :r]
            emit(f"rank_f{f}", *best_over_starts(
                lambda a, bb_, _B=B: A.fit(dhat, sd, i, j, a, bb_, wpair=w0, basis=_B)))

        # ---------------------------------------------------------------- robust arms
        for s in CAUCHY:
            emit(f"cauchy_s{s}", *best_over_starts(
                lambda a, bb_, _s=s: A.fit(dhat, sd, i, j, a, bb_, wpair=w0,
                                           kind="cauchy", s=_s)))
        for s in WELSCH:
            emit(f"welsch_s{s}", *best_over_starts(
                lambda a, bb_, _s=s: A.fit(dhat, sd, i, j, a, bb_, wpair=w0,
                                           kind="welsch", s=_s)))
        emit("gnc_cauchy", *best_over_starts(
            lambda a, b: RB.fit_gnc(dhat, sd, i, j, a, b, kind="cauchy",
                                    schedule=(8.0, 4.0, 2.0, 1.0), maxiter=150)))

        # ---------------------------------------------------------------- combinations
        wt_term = w0 * term_weight(n, i, j, 2, 0.25)
        wt_damp = w0 * damage_weight(Jsup, Gp, w0, 2.0)
        wt_both = w0 * term_weight(n, i, j, 2, 0.25) * damage_weight(Jsup, Gp, w0, 2.0)
        emit("cauchy1_x_term", *best_over_starts(
            lambda a, b: A.fit(dhat, sd, i, j, a, b, wpair=wt_term, kind="cauchy", s=1.0)))
        emit("cauchy1_x_damp", *best_over_starts(
            lambda a, b: A.fit(dhat, sd, i, j, a, b, wpair=wt_damp, kind="cauchy", s=1.0)))
        emit("term_x_damp", *best_over_starts(
            lambda a, b: A.fit(dhat, sd, i, j, a, b, wpair=wt_both)))

        # ---------------------------------------------------------------- ORACLE ceilings
        e = np.concatenate([A.wrap(sq_phi - nphi), A.wrap(sq_psi - npsi)])
        for f in KILLF:
            r = max(1, int(round(f * 2 * n)))
            P = V_sv[:, :r]
            e2 = e - P @ (P.T @ e)                     # delete the LOUD component of the error
            th = np.concatenate([nphi, npsi]) + e2
            emit(f"ORACLE_kill_loud_{f}", th[:n], th[n:])
            #: NORM-PRESERVING twin.  Deleting a component also SHRINKS |e|, so the arm above
            #: confounds "better direction" with "smaller error" -- exactly the confound this
            #: whole workstream exists to separate.  Rescaling to the original norm gives the
            #: honest ceiling: the SAME angular error budget, rotated out of the loud subspace.
            ne2 = float(np.linalg.norm(e2))
            if ne2 > 1e-12:
                th2 = np.concatenate([nphi, npsi]) + e2 * (np.linalg.norm(e) / ne2)
                emit(f"ORACLE_kill_loud_norm_{f}", th2[:n], th2[n:])
        oph, ops = sq_phi.copy(), sq_psi.copy()
        oph[:2] = nphi[:2]; oph[-2:] = nphi[-2:]
        ops[:2] = npsi[:2]; ops[-2:] = npsi[-2:]
        emit("ORACLE_term_native", oph, ops)

        row["align_start"] = float(A.alignment(
            A.sup_jacobian(nphi, npsi)[0],
            np.concatenate([A.wrap(np.asarray(S[0][0], float) - nphi),
                            A.wrap(np.asarray(S[0][1], float) - npsi)])))
        rows.append(row)

        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            with open(PATH, "w") as fh:
                json.dump({"rows": rows, "n_done": len(rows), "n_expected": len(tg)}, fh)
            el = time.time() - t_start
            print(f"  {c+1}/{len(tg)} checkpointed  ({el/60:.1f} min, "
                  f"{el/(c+1):.1f} s/target)", flush=True)

    return summarise(rows, tg)


# ------------------------------------------------------------------ leave-fold-out assembly
SWEEPS = {
    "TERMINAL_lfo": [f"term{b}_t{t}" for b, t in TERM],
    "JACW_lfo": [f"jacw_a{a}" for a in JACW],
    "DAMP_lfo": [f"damp_k{k}" for k in DAMP],
    "RANK_lfo": [f"rank_f{f}" for f in RANKF],
    "CAUCHY_lfo": [f"cauchy_s{s}" for s in CAUCHY],
    "ALL_ALIGN_lfo": ([f"term{b}_t{t}" for b, t in TERM] + [f"jacw_a{a}" for a in JACW] +
                      [f"damp_k{k}" for k in DAMP] + [f"rank_f{f}" for f in RANKF]),
}


def _lfo(rows, members, key="rmsd"):
    """Leave-fold-out arm: for held-out fold f use the member best on the OTHER folds."""
    folds = np.asarray([r["fold"] for r in rows], int)
    M = np.asarray([[r["arms"][m][key] for m in members] for r in rows], float)
    out = np.empty(len(rows))
    infold = np.empty(len(rows))
    pick = {}
    for f in sorted(set(folds.tolist())):
        te = folds == f
        tr = ~te
        k = int(M[tr].mean(0).argmin())
        out[te] = M[te, k]
        ki = int(M[te].mean(0).argmin())
        infold[te] = M[te, ki]
        pick[int(f)] = members[k]
    return out, infold, pick


def summarise(rows, tg=None):
    from s14 import ladder as LD
    pdbs = [r["pdb"] for r in rows]
    folds = np.asarray([r["fold"] for r in rows], int)
    fail = np.isin(pdbs, I.FAIL18)
    inc = LD.incumbent_rmsd()
    ref = np.asarray([inc[p] for p in pdbs], float)
    names = sorted(rows[0]["arms"].keys())
    keys = ("rmsd", "tors_rms_deg", "dist_mae", "resid_medz", "align", "align_self")

    arms = {}
    for a in names:
        arms[a] = {k: np.asarray([r["arms"][a][k] for r in rows], float) for k in keys}
    picks = {}
    for tag, mem in SWEEPS.items():
        mem = [m for m in mem if m in names]
        if not mem:
            continue
        v, infold, pk = _lfo(rows, mem)
        picks[tag] = pk
        arms[tag] = {"rmsd": v}
        arms[tag + "_INFOLD_ORACLE"] = {"rmsd": infold}
        for k in keys[1:]:
            #: carry the other axes from the fold-selected member
            M = {m: np.asarray([r["arms"][m][k] for r in rows], float) for m in mem}
            arms[tag][k] = np.asarray([M[pk[int(f)]][ii]
                                       for ii, f in enumerate(folds)], float)
            arms[tag + "_INFOLD_ORACLE"][k] = arms[tag][k]

    base = arms["squared"]["rmsd"]
    out = {"n": len(rows), "incumbent": float(ref.mean()), "picks": picks,
           "rows": rows, "arms": {}}
    for a in arms:
        v = arms[a]["rmsd"]
        out["arms"][a] = {
            **I.summary(v), "median": float(np.median(v)),
            "FAIL18": float(v[fail].mean()),
            "frac_under_2_5": float((v < 2.5).mean()),
            **{k: float(np.mean(arms[a][k])) for k in keys[1:] if k in arms[a]},
            "vs_squared": L.report(v, base, folds=folds, names=pdbs),
            "vs_incumbent": I.paired(v, ref, folds=folds, names=pdbs)}
    with open(PATH, "w") as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    I.write("s15_align_fit", out, n_expected=len(rows))
    report(out)
    return out


def report(out):
    print(f"\nn = {out['n']}   incumbent {out['incumbent']:.3f}\n")
    print(f"{'arm':<26}{'RMSD':>7}{'med':>7}{'<2.5':>6}{'FAIL18':>8}"
          f"{'tors':>7}{'dMAE':>7}{'medz':>7}{'align':>7}{'algS':>7}"
          f"{'   vs squared (RMSD)':<26}")
    order = sorted(out["arms"], key=lambda a: out["arms"][a]["mean"])
    for a in order:
        s = out["arms"][a]
        v = s["vs_squared"]
        print(f"{a:<26}{s['mean']:>7.3f}{s['median']:>7.3f}{s['frac_under_2_5']:>6.2f}"
              f"{s['FAIL18']:>8.3f}{s['tors_rms_deg']:>7.1f}{s['dist_mae']:>7.3f}"
              f"{s['resid_medz']:>7.3f}{s['align']:>7.3f}{s['align_self']:>7.3f}"
              f"   {v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]"
              f" W/L {v['n_better']}/{v['n_worse']}")
    print("\npicks (leave-fold-out):")
    for tag, pk in out["picks"].items():
        print(f"  {tag:<18} " + ", ".join(f"f{f}={m}" for f, m in sorted(pk.items())))


if __name__ == "__main__":
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else None
    run(limit=lim)
