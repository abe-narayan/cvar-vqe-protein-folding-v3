"""SPRINT 21 / WORKSTREAM B -- runner.

    python -m s21.qb3_run a  [n]     BLOCK A -- the ansatz zoo against its bond-dimension limit
    python -m s21.qb3_run q  [n]     BLOCK Q -- optimiser battery incl. Fisher natural gradient
    python -m s21.qb3_run e  [n]     BLOCK E -- the encoding as a GAUGE question (spsa, adam_fd)
    python -m s21.qb3_run en [n]     BLOCK E -- the same gauge sweep on Nelder-Mead
    python -m s21.qb3_run gate       the F-E5 harness gate: best_of_N bit-identity across charts

Hamiltonians, standardisation, starts, budget accounting and the RMSD readout are Sprint 20's,
imported unchanged.  Checkpointed per target; re-running resumes.  Every Ca-RMSD is ORACLE and
post-hoc, on the BUILT CHAIN basis, and is the argmin over everything evaluated.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

from s20 import qb2_lib as L
from s20 import qb2_opt as OP
from s20 import qb2_run as R20
from s21 import qb3_lib as K
from s15 import seed as SD

SEEDS = (0, 1, 2, 3)                 # 4 seeds, mandatory (s20 L14)
B = 512                              # evaluation budget per (target, objective, arm, seed)
N_T = 20

# ---- BLOCK A
A_KINDS = ("LEG", "AMBc", "DIST")
A_ANSATZ = ("prod", "mps_L1", "mps_L2", "mps_L3", "mps_L4", "mps_L2_nofinal",
            "share_L2", "sv_ring_L3")
A_SHOTS = 64
A_ALPHA = 0.25
SV_MAX_N = 14                        # `sv_ring_L3` is exact-statevector; capped and reported

# ---- BLOCK Q
Q_KINDS = ("LEG", "AMBc", "DIST")
Q_OPTS = ("adam", "sgd", "qng", "qng_diag", "spsa_ansatz")
Q_SHOTS = 32                         # 16 gradient steps at B=512; stated, not matched to s20's 64
Q_ALPHA = 0.25

# ---- BLOCK S (PREREG addendum s6): shots x iterations = B, an 8x sweep of optimisation
# effort with the encoding, chart, dimension, probe size and budget all held exactly fixed.
S_KINDS = ("LEG", "AMBc", "DIST")
S_SHOTS = (16, 32, 64, 128)          # -> 32, 16, 8, 4 gradient steps at B = 512
S_ALPHA = 0.25

# ---- BLOCK E
E_KINDS = ("LEG", "AMBc", "DIST")
E_ARMS = ("spsa", "adam_fd")
EN_ARMS = ("nelder",)
E_CHARTS = ("th", "th_nowrap", "th_s05", "th_s2", "emb_r1", "emb_r05", "emb_r2", "emb_norm")


def log(tag, msg):
    with open(os.path.join(K.RESULTS, f"qb3_{tag}.log"), "a") as fh:
        fh.write(msg + "\n")
    print(msg, flush=True)


def ck_path(tag):
    return os.path.join(K.RESULTS, f"qb3_{tag}.json")


def ck_load(tag):
    p = ck_path(tag)
    if os.path.exists(p):
        with open(p) as fh:
            return json.load(fh)
    return {}


def ck_save(tag, d):
    p = ck_path(tag)
    tmp = p + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(d, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    for _ in range(30):
        try:
            os.replace(tmp, p)
            return
        except PermissionError:
            time.sleep(0.4)


def finish(tag, want_keys, got):
    """A completion flag that requires the FULL configuration, not the subset run.

    `want_keys` is the complete (target, objective, arm, seed) product; the flag is written
    only if every one of them is present, and the shortfall is recorded otherwise.
    """
    miss = [k for k in want_keys if k not in got]
    meta = {"n_expected": len(want_keys), "n_present": len(want_keys) - len(miss),
            "missing_examples": miss[:12]}
    if miss:
        with open(os.path.join(K.RESULTS, f"qb3_{tag}__PARTIAL_"), "w") as fh:
            json.dump(meta, fh)
        log(tag, f"[PARTIAL] {len(miss)}/{len(want_keys)} cells missing")
        return False
    with open(os.path.join(K.RESULTS, f"qb3_{tag}_COMPLETE"), "w") as fh:
        fh.write(time.strftime("%Y-%m-%d %H:%M:%S") + " " + json.dumps(meta))
    pp = os.path.join(K.RESULTS, f"qb3_{tag}__PARTIAL_")
    if os.path.exists(pp):
        os.remove(pp)
    log(tag, f"[COMPLETE] all {len(want_keys)} cells present")
    return True


def with_mem_retry(tag, fn, tries=4):
    """Run `fn`, surviving `core.amber.memory_guard`'s refusal to start an OpenMM context.

    The guard raises above 92% to protect the other processes on this shared box; raising the
    ceiling is not an option, so a refusal is waited out and retried rather than allowed to kill
    the whole block.  The number of refusals actually caught is returned and recorded -- a guard
    that never fires is not evidence (BRIEF section 7 rule 3).
    """
    caught = 0
    for k in range(tries):
        try:
            return fn(), caught
        except MemoryError as exc:
            caught += 1
            log(tag, f"[mem] guard refused ({exc.args[0][:60]}...), retry {k + 1}/{tries}")
            mem_wait(tag, limit=86.0, patience=900)
    log(tag, "[mem] guard refused on every retry; target skipped and recorded")
    return None, caught


def mem_wait(tag, limit=90.5, patience=600):
    """Wait until the box is below `limit` percent before starting OpenMM contexts.

    `core.amber.memory_guard` refuses above 92% to protect the other processes on this shared
    box.  Raising that ceiling is not an option, so the lane waits.  Every wait is logged, and
    the number of times the wait actually FIRED is reported -- a guard that never fires is not
    evidence (BRIEF section 7 rule 3).

    CORRECTED MID-RUN, and recorded: the first setting was `limit=88, patience=3600`, which is
    3.5 points BELOW the ceiling that actually refuses.  On a box hovering at 89% that spends up
    to an hour of wall clock per target waiting for a guard that was never going to fire -- a
    conservative threshold that costs more than the thing it protects against.  The wait is now
    set just under the real refusal point, and `with_mem_retry` handles an actual refusal.
    """
    from core import amber as am
    waited = 0
    fired = 0
    while am.memory_percent() > limit and waited < patience:
        if fired == 0:
            log(tag, f"[mem] {am.memory_percent():.0f}% > {limit:.0f}%, waiting")
        fired += 1
        time.sleep(20)
        waited += 20
    return fired


def _field(kind, tgt, budget, sp, cal_k, keep=False):
    h = L.Ham(kind, tgt, budget=budget, keep=keep, sp=sp)
    h.mu_ref = cal_k["mu_ref"]; h.sd_ref = cal_k["sd_ref"]
    h.med_ref = cal_k["med_ref"]; h.iqr_ref = cal_k["iqr_ref"]
    F = OP.Field(h, tgt)
    F.set_pen(cal_k.get("ref_std_max", 6.0))
    return F, h


def _readout(F, h, tgt, z0):
    """The single, shared readout: argmin over everything evaluated, on the BUILT CHAIN."""
    rm, braw, disp = float("nan"), float("nan"), float("nan")
    if F.best_z is not None:
        rm = float(L.rmsd_of(F.best_z[None], tgt)[0])
        braw = float(h.raw(*L.unpack(F.best_z[None], tgt["n"]))[0])
        disp = float(np.linalg.norm(L.wrap(F.best_z - np.asarray(z0, float))))
    return {"best_std": float(F.best), "best_raw": braw, "rmsd_ORACLE": rm,
            "disp_ang": disp, "used": int(h.used), "n_nonfinite": int(F.n_nonfinite)}


XREAD = ("DIST", "LEG")          # readout Hamiltonians that cost no OpenMM call


def cross_readout(h, tgt, sp=None):
    """BRIEF section 2, open question (b): a READOUT Hamiltonian different from the TRAINING one.

    CVaR is the training objective; the readout is an argmin (the BRIEF's own mid-sprint
    correction).  Nothing forces the argmin to be taken under the SAME H the sampler trained on.
    Every arm already evaluated a set of configurations, so re-taking the argmin under a
    different H is FREE -- no extra budget, no extra structure, and for `DIST`/`LEG` no OpenMM
    call at all (measured: 0.24 s for both over 512 configurations).

    Returns the Ca-RMSD of each readout's argmin over the SAME evaluated set, so the training
    Hamiltonian and the readout Hamiltonian are separated experimentally rather than assumed
    equal.  `xread_<K>` where K == the training kind is an IDENTITY CHECK and must reproduce the
    arm's own `rmsd_ORACLE`.
    """
    phi, psi, _e = h.seen()
    out = {}
    if len(phi) == 0:
        return out
    for rk in XREAD:
        hr = L.Ham(rk, tgt, budget=10 ** 9, keep=False, sp=sp)
        er = np.asarray(hr.raw(phi, psi), float)
        ok = np.isfinite(er)
        if not ok.any():
            continue
        k = int(np.flatnonzero(ok)[int(np.argmin(er[ok]))])
        z = L.pack(phi[k:k + 1], psi[k:k + 1])
        out[f"xread_{rk}_rmsd"] = float(L.rmsd_of(z, tgt)[0])
    out["xread_n_seen"] = int(len(phi))
    return out


# ============================================================ BLOCK A / Q  (variational)
def run_vqe_block(pdbs, tag, kinds, rows_spec, budget=B, full_kinds=None):
    """`rows_spec(tgt)` -> list of (row_key, kwargs for `arm_vqe2`) or ("best_of_N", None).

    `kinds` is what THIS call runs (objectives can be run in separate passes when the box's
    memory guard blocks OpenMM); `full_kinds` is the FULL declared configuration the completion
    flag requires.  A flag is never written for the subset a call happened to run.
    """
    full_kinds = tuple(full_kinds or kinds)
    out = ck_load(tag)
    for ti, pdb in enumerate(pdbs):
        rec = out.get(pdb)
        if rec is not None and all(
                "xread_n_seen" in rec["rows"].get(f"{k}|{key}|{s}", {})
                for k in kinds for key, _ in rows_spec(rec) for s in SEEDS):
            continue
        if any(k in ("AMB", "AMBc") for k in kinds):
            R20.mem_hold(tag=tag)
            mem_wait(tag)
        t0 = time.time()
        tgt = L.target(pdb)
        Z = R20.starts_for(tgt)
        got, n_refused = with_mem_retry(tag, lambda: R20.prep(tgt, kinds))
        if got is None:
            continue
        H, sp, cal, ver = got
        if rec is None:
            rec = {"n": tgt["n"], "fold": tgt["fold"], "budget": budget, "rows": {}}
        rec["mem_guard_refusals"] = rec.get("mem_guard_refusals", 0) + n_refused
        rec["amber_exact_maxrel"] = ver[0] if ver else rec.get("amber_exact_maxrel")
        rec["amber_exact_ncmp"] = ver[1] if ver else rec.get("amber_exact_ncmp")
        for k in kinds:
            for key, kw in rows_spec(tgt):
                for s in SEEDS:
                    if f"{k}|{key}|{s}" in rec["rows"]:
                        continue
                    z0 = Z[s % len(Z)]
                    F, h = _field(k, tgt, budget, sp, cal[k], keep=True)
                    rng = SD.stable_rng(pdb, k, key, s, salt=K.SALT)
                    tA = time.time()
                    if kw is None:
                        info = OP.arm_bestofn(F, z0, rng)
                    else:
                        info = K.arm_vqe2(F, z0, rng, **kw)
                    row = {**_readout(F, h, tgt, z0), "wall": time.time() - tA, **info}
                    row.update(cross_readout(h, tgt, sp))
                    rec["rows"][f"{k}|{key}|{s}"] = row
        rec["ref"] = {
            "start_rmsd_ORACLE": float(L.rmsd_of(Z[0][None], tgt)[0]),
            "pool_best_ORACLE": float(np.min(L.rmsd_of(L.pack(tgt["PHI"], tgt["PSI"]), tgt))),
            "pool_mean_ORACLE": float(np.mean(L.rmsd_of(L.pack(tgt["PHI"], tgt["PSI"]), tgt))),
        }
        out[pdb] = rec
        ck_save(tag, out)
        best = {k: np.nanmin([v["rmsd_ORACLE"] for kk, v in rec["rows"].items()
                              if kk.startswith(k + "|")]) for k in kinds}
        log(tag, f"[{ti+1}/{len(pdbs)}] {pdb} n={tgt['n']} {time.time()-t0:.0f}s "
                 + " ".join(f"{k}:{v:.2f}" for k, v in best.items()))
    # FULL-configuration flag: every target x objective x arm x seed the spec calls for at
    # that target's own n (the `sv_ring_L3` cap is part of the declared configuration).
    want = []
    for pdb in pdbs:
        if pdb not in out:
            want.append((pdb, "<target missing>"))
            continue
        for k in full_kinds:
            for key, _kw in rows_spec({"n": out[pdb]["n"], "pdb": pdb}):
                for s in SEEDS:
                    want.append((pdb, f"{k}|{key}|{s}"))
    got = {(p, r) for p, v in out.items() for r in v.get("rows", {})}
    finish(tag, want, got)
    return out


def spec_A(tgt):
    rows = []
    for a in A_ANSATZ:
        if a == "sv_ring_L3" and int(tgt.get("n", 99)) > SV_MAX_N:
            continue
        rows.append((a, {"ansatz": a, "opt": "adam", "alpha": A_ALPHA, "shots": A_SHOTS,
                         "train": True}))
    rows.append(("untrained_L2", {"ansatz": "mps_L2", "opt": "adam", "alpha": A_ALPHA,
                                  "shots": A_SHOTS, "train": False}))
    rows.append(("best_of_N", None))
    return rows


def spec_S(tgt):
    rows = [(f"sh{sh}", {"ansatz": "mps_L2", "opt": "adam", "alpha": S_ALPHA, "shots": sh,
                         "train": True}) for sh in S_SHOTS]
    rows.append(("best_of_N", None))
    return rows


def spec_Q(tgt):
    rows = [(f"opt_{o}", {"ansatz": "mps_L2", "opt": o, "alpha": Q_ALPHA, "shots": Q_SHOTS,
                          "train": True}) for o in Q_OPTS]
    rows.append(("untrained_L2", {"ansatz": "mps_L2", "opt": "adam", "alpha": Q_ALPHA,
                                  "shots": Q_SHOTS, "train": False}))
    rows.append(("best_of_N", None))
    return rows


# ============================================================ BLOCK E  (the gauge sweep)
def run_enc(pdbs, tag, arms, kinds=E_KINDS, charts=E_CHARTS, budget=B, full_kinds=None):
    full_kinds = tuple(full_kinds or kinds)
    out = ck_load(tag)
    want = [(p, f"{k}|{a}|{c}|{s}") for p in pdbs for k in full_kinds for a in arms
            for c in charts for s in SEEDS]
    for ti, pdb in enumerate(pdbs):
        rec = out.get(pdb)
        if rec is not None and all(f"{k}|{a}|{c}|{s}" in rec["rows"] for k in kinds
                                   for a in arms for c in charts for s in SEEDS):
            continue
        if any(k in ("AMB", "AMBc") for k in kinds):
            R20.mem_hold(tag=tag)
            mem_wait(tag)
        t0 = time.time()
        tgt = L.target(pdb)
        Z = R20.starts_for(tgt)
        got, n_refused = with_mem_retry(tag, lambda: R20.prep(tgt, kinds))
        if got is None:
            continue
        H, sp, cal, ver = got
        if rec is None:
            rec = {"n": tgt["n"], "fold": tgt["fold"], "budget": budget, "rows": {}}
        rec["mem_guard_refusals"] = rec.get("mem_guard_refusals", 0) + n_refused
        rec["amber_exact_maxrel"] = ver[0] if ver else rec.get("amber_exact_maxrel")
        for k in kinds:
            for arm in arms:
                for cname in charts:
                    for s in SEEDS:
                        if f"{k}|{arm}|{cname}|{s}" in rec["rows"]:
                            continue
                        z0 = Z[s % len(Z)]
                        F, h = _field(k, tgt, budget, sp, cal[k])
                        C = K.CHARTS[cname](F)
                        x0 = C.from_z(z0[None])[0]
                        rng = SD.stable_rng(pdb, k, arm, cname, s, salt=K.SALT)
                        tA = time.time()
                        # `track_quality` spends 3 central-FD probes of 2*d units each, and
                        # `d` DOUBLES in the circle chart -- 288 of 512 budget units embedded
                        # against 144 in the angle chart.  The audit lane showed that single
                        # diagnostic accounted for the WHOLE of SPSA's 159-vs-208 iteration
                        # deficit in Sprint 20.  It is switched off here so the gauge sweep
                        # is matched in the space the operator works in (BRIEF s7 rule 1).
                        info = (OP.arm_spsa(C, x0, rng, track_quality=False)
                                if arm == "spsa" else OP.ARMS[arm](C, x0, rng))
                        rec["rows"][f"{k}|{arm}|{cname}|{s}"] = {
                            **_readout(F, h, tgt, z0), "wall": time.time() - tA,
                            **C.info(), **{kk: vv for kk, vv in info.items()
                                           if not isinstance(vv, list)}}
        out[pdb] = rec
        ck_save(tag, out)
        log(tag, f"[{ti+1}/{len(pdbs)}] {pdb} n={tgt['n']} {time.time()-t0:.0f}s")
    got = {(p, r) for p, v in out.items() for r in v.get("rows", {})}
    finish(tag, want, got)
    return out


# ============================================================ the F-E5 harness gate
def run_gate(pdbs):
    """`best_of_N` must emit a BIT-IDENTICAL torsion set in every chart at the same seed.

    It draws from the same von Mises basin mixtures and is mapped back through the same
    retraction, so any difference is a harness bug, not physics.  The gate reports the measured
    max |delta| and the number of comparisons -- a gate that never fires is not evidence.
    """
    rows = {}
    for pdb in pdbs:
        tgt = L.target(pdb)
        base = None
        worst = 0.0
        ncmp = 0
        for cname in E_CHARTS:
            h = L.Ham("DIST", tgt, budget=256, keep=True)
            h.calibrate(tgt["PHI"], tgt["PSI"])
            F = OP.Field(h, tgt)
            F.set_pen(6.0)
            C = K.CHARTS[cname](F)
            rng = SD.stable_rng(pdb, "gate", salt=K.SALT)
            OP.arm_bestofn(C, C.from_z(R20.starts_for(tgt)[0][None])[0], rng)
            phi, psi, e = h.seen()
            zz = L.pack(phi, psi)
            if base is None:
                base = zz
            else:
                m = min(len(base), len(zz))
                worst = max(worst, float(np.max(np.abs(base[:m] - zz[:m]))))
                ncmp += m
        rows[pdb] = {"max_abs_delta": worst, "n_compared": ncmp,
                     "n_charts": len(E_CHARTS)}
        log("gate", f"{pdb} max|delta|={worst:.3e} over {ncmp} configs")
    # The flag must require the FULL configuration, not the subset this call ran.
    full = len(rows) == N_T
    K.write21("qb3_gate", {"charts": list(E_CHARTS), "rows": rows,
                           "n_targets": len(rows), "n_targets_required": N_T,
                           "worst_overall": max(v["max_abs_delta"] for v in rows.values()),
                           "total_compared": sum(v["n_compared"] for v in rows.values())},
              complete=full)
    if not full:
        with open(os.path.join(K.RESULTS, "qb3_gate__PARTIAL_"), "w") as fh:
            json.dump({"n_targets": len(rows), "n_targets_required": N_T}, fh)
    return rows


# ============================================================ main
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "a"
    n_t = int(sys.argv[2]) if len(sys.argv) > 2 else N_T
    ss = L.subset(N_T)[:n_t]
    pdbs = [t["pdb"] for t in ss]
    K.write21("qb3_config", {"salt": K.SALT, "subset": pdbs, "B": B, "SEEDS": list(SEEDS),
                             "A_ANSATZ": list(A_ANSATZ), "A_SHOTS": A_SHOTS,
                             "A_ALPHA": A_ALPHA, "SV_MAX_N": SV_MAX_N,
                             "Q_OPTS": list(Q_OPTS), "Q_SHOTS": Q_SHOTS, "Q_ALPHA": Q_ALPHA,
                             "E_CHARTS": list(E_CHARTS), "E_ARMS": list(E_ARMS),
                             "kinds_A": list(A_KINDS), "kinds_Q": list(Q_KINDS),
                             "S_SHOTS": list(S_SHOTS), "S_ALPHA": S_ALPHA,
                             "kinds_S": list(S_KINDS),
                             "kinds_E": list(E_KINDS), "mode": mode,
                             "n_targets": len(pdbs)}, complete=True)
    # `mode` may carry an objective suffix -- `a:DIST` runs only that objective's rows.  The
    # completion flag still requires the FULL declared kind list either way.
    sub = None
    if ":" in mode:
        mode, sub = mode.split(":", 1)
    def _k(full):
        return (sub,) if sub else full
    if mode == "a":
        run_vqe_block(pdbs, "a", _k(A_KINDS), spec_A, full_kinds=A_KINDS)
    elif mode == "q":
        run_vqe_block(pdbs, "q", _k(Q_KINDS), spec_Q, full_kinds=Q_KINDS)
    elif mode == "s":
        run_vqe_block(pdbs, "s", _k(S_KINDS), spec_S, full_kinds=S_KINDS)
    elif mode == "e":
        run_enc(pdbs, "e", E_ARMS, kinds=_k(E_KINDS), full_kinds=E_KINDS)
    elif mode == "en":
        run_enc(pdbs, "en", EN_ARMS, kinds=_k(E_KINDS), full_kinds=E_KINDS)
    elif mode == "gate":
        run_gate(pdbs)
    else:
        raise SystemExit(f"unknown mode {mode}")
    print("done")


if __name__ == "__main__":
    main()
