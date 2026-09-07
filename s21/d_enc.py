"""s21/d_enc.py -- D1: IS THE ENCODING LEVER A STEP-COUNT EFFECT?

    python -m s21.d_enc gate      the bit-for-bit soundness gate against Sprint 20's own rows
    python -m s21.d_enc run       the budget grid
    python -m s21.d_enc report    the analysis

WHAT IS UNDER ATTACK.  s21/BRIEF.md section 5 item 3 and Sprint 20 section 8: `theta` versus
`(cos theta, sin theta)` -- physically identical, because `arctan2` is scale-invariant so the
radius is pure gauge -- won on final Ca-RMSD in 11 of 12 cells, p = 0.006, up to -0.649 A.

WHAT THE AUDIT FOUND BEFORE THIS FILE WAS WRITTEN (in `s20/results/qb2_enc.json` and
`qb2_opt.json`, no rerun needed, reported as audit rather than as experiment):

  * `AMB|nelder` and `AMBc|nelder` are BIT-IDENTICAL on all 10 targets x 2 seeds in BOTH
    encodings.  `AMBc = sign(E)*log1p(|E|)` is strictly monotone and Nelder-Mead is
    comparison-only, so this is an IDENTITY.  Two of the "twelve cells" are one measurement.
  * At the brief's mandated unit (target), the effect nonetheless survives:
    -0.2167 [-0.3239, -0.1095], 9W/1L over 10 targets.  The MEASUREMENT is real.
  * THE EMBEDDING TAKES FEWER OPTIMISATION STEPS IN EVERY CELL.  From the artefacts' own
    `iters`: adam_fd 5.1 vs 10.5 (ratio 0.49), spsa 159.2 vs 207.6 (0.77).  And spsa's deficit
    is not the declared 4n-vs-2n dimension effect -- it is `arm_spsa`'s `track_quality`
    DIAGNOSTIC, three central-FD probes of 2*d units each, and `d` doubles in the embedding, so
    the diagnostic eats 288 of 512 budget units embedded against 144 in theta.

WHY THAT MATTERS HERE SPECIFICALLY.  On this instrument, optimising the deployed objectives
HARDER makes the STRUCTURE worse -- Q3 (the arm optimising hardest is second-worst on the
board), `search-saturates-discrimination-binds`, `concentration-is-wrong-when-discrimination-
binds`.  The start is a retrieval-derived good structure.  So "took fewer steps" and "found a
better structure" are the same event here, and the Sprint-20 comparison cannot separate them.
Its own objective column agrees: the embedding reaches a WORSE objective in 9 of 12 cells.

THE CONTROL THAT DECIDES IT.  Matched realised STEP COUNT, not matched budget -- BRIEF section 7
rule 1, a control must be matched in the space the operator works in.  Rather than pick one
handicapping budget, the step-count response curve is MEASURED on a grid in both encodings and
the matched comparison is read off it:

    encodings   theta | emb            budgets   256 | 512 | 1024
    kinds       LEG (free) | AMBc (the biggest Sprint-20 effect, -0.649)
    arms        spsa | adam_fd         seeds     0 1 2 3   (BRIEF section 6 minimum; s20 used 2)

`track_quality` is OFF in every arm here, which removes the asymmetric diagnostic cost, so the
two mechanisms are separated: if the gap dies with the diagnostic off, the diagnostic was the
whole thing; if it survives that but dies under step-matching, step count is.

TWO INSTRUMENTS THE SPRINT-20 ARTEFACTS DO NOT CARRY.

    disp_best / disp_final   ||wrap(z - z0)|| in TORSION space, radians -- how far the arm
                             actually moved.  The Sprint-20 rows record `update_norm` in
                             whichever space the arm ran in, which is a basis mismatch between
                             the two arms of the headline comparison.
    r_first / r_last         the gauge radius ||u_i|| per angle, embedded only.  `arctan2` is
                             scale-invariant so the objective's gradient has EXACTLY zero radial
                             component; a finite-difference optimiser injects noise into it, and
                             `adam_fd` normalises per coordinate by sqrt(v) so it takes a
                             FULL-SIZE step along a direction whose gradient is pure noise.  If
                             r drifts up, the effective angular step decays as 1/r -- an
                             unintended annealing schedule, which is mechanism (b).

PRE-REGISTERED IN `s21/PREREG_D.md` D1, INCLUDING THE FALSIFIER OF MY OWN HYPOTHESIS: if at
matched step count the embedding still wins by more than the MDE with a CI excluding zero in
BOTH matching directions, the lever is not a step-count effect and I withdraw D1.

SOUNDNESS GATE, run before any number is read.  `arm_spsa_rec` / `arm_adam_rec` below are
VERBATIM copies of `s20.qb2_opt.arm_spsa` / `arm_adam` with recording added and nothing else
changed.  With `track_quality=True`, the same seed, start and budget they must reproduce Sprint
20's stored `rmsd_ORACLE` and `best_std` BIT-FOR-BIT.  The gate's firing count is reported --
a gate that never fires is not evidence (BRIEF section 7 rule 3).
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s20 import qb2_lib as L          # noqa: E402
from s20 import qb2_opt as OP         # noqa: E402
from s20 import qb2_run as RUN        # noqa: E402
from s15 import seed as SD            # noqa: E402
from core import quantum as Q         # noqa: E402

#: the SAME 10 targets Sprint 20's encoding panel used, so the attack is paired to its own
#: instrument rather than to a fresh sample.
TARGETS = ["1CS9", "1ID6", "1TOR", "2JN5", "2LU6", "2MFV", "2MIG", "2MK7", "2MSA", "2O0S"]
KINDS = ("LEG", "AMBc")
ARMS = ("spsa", "adam_fd")
ENCS = ("theta", "emb")
BUDGETS = (256, 512, 1024)
SEEDS = (0, 1, 2, 3)
CFG = {"targets": TARGETS, "kinds": list(KINDS), "arms": list(ARMS), "encs": list(ENCS),
       "budgets": list(BUDGETS), "seeds": list(SEEDS), "track_quality": False}
N_CELLS = len(TARGETS) * len(KINDS) * len(ARMS) * len(ENCS) * len(BUDGETS) * len(SEEDS)
GATE_TARGETS = TARGETS[:3]
TAG = "d_enc"          # set from argv so a LEG-only pass and an AMBc pass never share a file


def _prep_kinds():
    """Only build the OpenMM context when an AMBER kind is actually in the configuration --
    the box is shared and `builder_for` refuses above a 92% memory ceiling."""
    k = list(KINDS)
    if "AMBc" in k and "AMB" not in k:
        k = ["AMB"] + k          # AMBc reuses AMB's pool single points; prep() shares them
    return tuple(k)


def _mem_hold(min_gb=1.6, max_wait=5400):
    """The box is shared and `core.amber.builder_for` refuses above a 92% ceiling.  Wait rather
    than crash, and never push the machine past the brief's limit."""
    from s12 import instrument as I
    if not any(k.startswith("AMB") for k in KINDS):
        return 0
    from core import amber as am
    waited = 0
    while (I.free_gb() < min_gb or am.memory_percent() > 88.0) and waited < max_wait:
        print(f"  [wait] free {I.free_gb():.2f} GB, mem {am.memory_percent():.0f}%", flush=True)
        time.sleep(30); waited += 30
    return waited


# ---------------------------------------------------------------- instrumented arms
def arm_spsa_rec(F, z0, rng, a=0.20, c=L.FD_H, A_frac=0.1, track_quality=False, to_z=None):
    """VERBATIM `s20.qb2_opt.arm_spsa`, with recording added and nothing else changed."""
    z = np.asarray(z0, float).copy()
    it = 0
    max_it = max(1, F.left // 2)
    A = max(1.0, A_frac * max_it)
    cos, unorm = [], []
    checks = {int(max_it * f) for f in (0.05, 0.4, 0.8)}
    while F.left >= 2:
        it += 1
        ak = a / (it + A) ** 0.602
        ck = c / it ** 0.101
        g, _d = OP._grad_spsa(F, z, ck, rng)
        if g is None:
            break
        if track_quality and it in checks and F.left >= 4 * F.n + 2:
            gt = OP._grad_fd(F, z)
            if gt is not None and np.linalg.norm(gt) > 0 and np.linalg.norm(g) > 0:
                cos.append(float(g @ gt / (np.linalg.norm(g) * np.linalg.norm(gt))))
        step = ak * g
        unorm.append(float(np.linalg.norm(step)))
        z = F.wrap(z - step)
    return {"iters": it, "spsa_cos_fd": cos, "z_final": z,
            "update_norm_med": float(np.median(unorm)) if unorm else float("nan"),
            "update_norm_p95": float(np.percentile(unorm, 95)) if unorm else float("nan")}


def arm_adam_rec(F, z0, rng, lr=0.10, to_z=None):
    """VERBATIM `s20.qb2_opt.arm_adam`, with recording added and nothing else changed."""
    z = np.asarray(z0, float).copy()
    opt = Q.Adam(len(z), lr=lr)
    it, unorm = 0, []
    while F.left >= 4 * F.n:
        g = OP._grad_fd(F, z)
        if g is None:
            break
        it += 1
        zn = opt.step(z, g)
        unorm.append(float(np.linalg.norm(zn - z)))
        z = F.wrap(zn)
    return {"iters": it, "z_final": z,
            "update_norm_med": float(np.median(unorm)) if unorm else float("nan"),
            "update_norm_p95": float(np.percentile(unorm, 95)) if unorm else float("nan")}


RECARMS = {"spsa": arm_spsa_rec, "adam_fd": arm_adam_rec}


def _radius(u, nz2):
    """Gauge radius per angle for an embedded coordinate vector: r_i = hypot(cos_i, sin_i)."""
    u = np.asarray(u, float)
    return np.hypot(u[:nz2], u[nz2:])


def run_cell(kind, tgt, arm, seed, z0, budget, sp, cal_k, enc, track_quality=False):
    """One (kind, arm, seed, budget, encoding) cell.  Mirrors `s20.qb2_run._run_arm`."""
    h = L.Ham(kind, tgt, budget=budget, keep=False, sp=sp)
    h.mu_ref = cal_k["mu_ref"]; h.sd_ref = cal_k["sd_ref"]
    h.med_ref = cal_k["med_ref"]; h.iqr_ref = cal_k["iqr_ref"]
    F = OP.Field(h, tgt)
    F.set_pen(cal_k.get("ref_std_max", 6.0))
    #: IDENTICAL rng key to Sprint 20's, so the theta@512 cell reproduces its row exactly
    rng = SD.stable_rng(tgt["pdb"], kind, arm, seed, salt=L.SALT)
    t0 = time.time()
    rec = {}
    if enc == "emb":
        E = RUN._EmbField(F)
        nz2 = 2 * E.nz                       # number of ANGLES (phi and psi per residue)
        u0 = np.concatenate([np.cos(z0), np.sin(z0)])
        info = RECARMS[arm](E, u0, rng, track_quality=track_quality) if arm == "spsa" \
            else RECARMS[arm](E, u0, rng)
        u_end = np.asarray(info.pop("z_final"), float)
        z_end = E.to_z(u_end[None])[0]
        rec["r_first"] = float(np.mean(_radius(u0, nz2)))
        rec["r_last"] = float(np.mean(_radius(u_end, nz2)))
        rec["r_last_max"] = float(np.max(_radius(u_end, nz2)))
    else:
        info = RECARMS[arm](F, z0, rng, track_quality=track_quality) if arm == "spsa" \
            else RECARMS[arm](F, z0, rng)
        z_end = np.asarray(info.pop("z_final"), float)
        rec["r_first"] = rec["r_last"] = rec["r_last_max"] = 1.0
    rm, braw = float("nan"), float("nan")
    disp_best = float("nan")
    if F.best_z is not None:
        rm = float(L.rmsd_of(F.best_z[None], tgt)[0])
        braw = float(h.raw(*L.unpack(F.best_z[None], tgt["n"]))[0])
        disp_best = float(np.linalg.norm(L.wrap(F.best_z - z0)))
    rec.update({"best_std": float(F.best), "best_raw": braw, "rmsd_ORACLE": rm,
                "used": int(h.used), "n_nonfinite": int(F.n_nonfinite),
                "wall": time.time() - t0,
                "disp_best": disp_best,
                "disp_final": float(np.linalg.norm(L.wrap(z_end - z0))),
                **{k: v for k, v in info.items() if k != "spsa_cos_fd"}})
    return rec


# ---------------------------------------------------------------- soundness gate
def gate():
    """Reproduce Sprint 20's OWN stored rows bit-for-bit with the copied arms."""
    ref = json.load(open(os.path.join(ROOT, "s20", "results", "qb2_enc.json")))
    refo = json.load(open(os.path.join(ROOT, "s20", "results", "qb2_opt.json")))
    fired, bad = 0, []
    out = []
    for pdb in GATE_TARGETS:
        _mem_hold()
        tgt = L.target(pdb)
        Z = RUN.starts_for(tgt)
        H, sp, cal, ver = RUN.prep(tgt, _prep_kinds())
        for kind in KINDS:
            for arm in ARMS:
                for sd_ in (0, 1):
                    for enc, src in (("theta", refo), ("emb", ref)):
                        want = src[pdb]["rows"][f"{kind}|{arm}|{sd_}"]
                        got = run_cell(kind, tgt, arm, sd_, Z[sd_ % len(Z)], 512, sp,
                                       cal[kind], enc, track_quality=True)
                        fired += 1
                        ok = (got["rmsd_ORACLE"] == want["rmsd_ORACLE"]
                              and got["best_std"] == want["best_std"]
                              and got["iters"] == want["iters"])
                        if not ok:
                            bad.append({"pdb": pdb, "kind": kind, "arm": arm, "seed": sd_,
                                        "enc": enc,
                                        "got": [got["rmsd_ORACLE"], got["best_std"], got["iters"]],
                                        "want": [want["rmsd_ORACLE"], want["best_std"],
                                                 want["iters"]]})
                        out.append({"pdb": pdb, "kind": kind, "arm": arm, "seed": sd_,
                                    "enc": enc, "ok": bool(ok)})
        print(f"  gate {pdb}: {fired} comparisons, {len(bad)} mismatches", flush=True)
        json.dump({"n_comparisons": fired, "n_mismatch": len(bad), "mismatches": bad[:20],
                   "gate_targets": GATE_TARGETS, "rows": out},
                  open(os.path.join(RESULTS, "d_enc_gate.json"), "w"), indent=1)
    res = {"n_comparisons": fired, "n_mismatch": len(bad), "mismatches": bad[:20],
           "gate_targets": GATE_TARGETS,
           "PASS": len(bad) == 0 and fired == len(GATE_TARGETS) * len(KINDS) * len(ARMS)
           * len(ENCS) * 2, "rows": out}
    json.dump(res, open(os.path.join(RESULTS, "d_enc_gate.json"), "w"), indent=1)
    print(f"\nGATE: {fired} comparisons, {len(bad)} mismatches -> "
          f"{'PASS' if res['PASS'] else 'FAIL'}")
    if res["PASS"]:
        with open(os.path.join(RESULTS, "d_enc_gate.COMPLETE"), "w") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} n_comparisons={fired} mismatches=0\n")
    return res


# ---------------------------------------------------------------- the grid
def _ck():
    p = os.path.join(RESULTS, f"{TAG}.json")
    return json.load(open(p)) if os.path.exists(p) else {}


def _save(d):
    p = os.path.join(RESULTS, f"{TAG}.json")
    tmp = p + f".tmp{os.getpid()}"
    json.dump(d, open(tmp, "w"), default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, p)


def run():
    out = _ck()
    t00 = time.time()
    for ti, pdb in enumerate(TARGETS):
        if pdb in out and len(out[pdb]["rows"]) == len(KINDS) * len(ARMS) * len(ENCS) \
                * len(BUDGETS) * len(SEEDS):
            continue
        _mem_hold()
        t0 = time.time()
        tgt = L.target(pdb)
        Z = RUN.starts_for(tgt)
        H, sp, cal, ver = RUN.prep(tgt, _prep_kinds())
        rec = {"n": tgt["n"], "fold": tgt["fold"], "cal": cal,
               "amber_exact_maxrel": ver[0] if ver else None,
               "cfg": CFG, "rows": {}}
        for kind in KINDS:
            for arm in ARMS:
                for enc in ENCS:
                    for B in BUDGETS:
                        for sd_ in SEEDS:
                            r = run_cell(kind, tgt, arm, sd_, Z[sd_ % len(Z)], B, sp,
                                         cal[kind], enc, track_quality=False)
                            rec["rows"][f"{kind}|{arm}|{enc}|{B}|{sd_}"] = r
        out[pdb] = rec
        _save(out)
        print(f"[{ti+1}/{len(TARGETS)}] {pdb} n={tgt['n']} {time.time()-t0:.0f}s "
              f"(total {time.time()-t00:.0f}s)", flush=True)
    _finish(out)
    return out


def _finish(out):
    """COMPLETION FLAG.  Requires the FULL configuration -- every target, kind, arm, encoding,
    budget and seed present with a finite RMSD -- not the subset the run was called with.
    (s21/results/tailprice.json was stamped complete:true with 3/3 rows SKIPPED.)"""
    need = [f"{k}|{a}|{e}|{B}|{s}" for k in KINDS for a in ARMS for e in ENCS
            for B in BUDGETS for s in SEEDS]
    miss = []
    for p in TARGETS:
        if p not in out:
            miss.append((p, "ABSENT")); continue
        for key in need:
            r = out[p]["rows"].get(key)
            if r is None:
                miss.append((p, key))
            elif not np.isfinite(r.get("rmsd_ORACLE", np.nan)):
                miss.append((p, key + " NONFINITE"))
    if miss:
        print(f"INCOMPLETE: {len(miss)} missing/nonfinite cells, e.g. {miss[:5]}")
        return False
    with open(os.path.join(RESULTS, f"{TAG}.COMPLETE"), "w") as fh:
        fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                 f"targets={len(TARGETS)} kinds={list(KINDS)} arms={list(ARMS)} "
                 f"encs={list(ENCS)} budgets={list(BUDGETS)} seeds={list(SEEDS)} "
                 f"cells={N_CELLS} track_quality=False\n")
    print(f"COMPLETE: all {N_CELLS} cells present and finite.")
    return True


# ---------------------------------------------------------------- report
def _boot(d, rng, B=20000):
    d = np.asarray(d, float); k = len(d)
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def _signp(n, k):
    s = sum(math.comb(n, i) for i in range(0, min(k, n - k) + 1))
    return min(1.0, 2 * s / 2 ** n)


def report():
    d = _ck()
    P = [p for p in TARGETS if p in d]
    rng = np.random.default_rng(20210907)
    folds = np.array([d[p]["fold"] for p in P])

    def cell(kind, arm, enc, B, field="rmsd_ORACLE"):
        return np.array([np.mean([d[p]["rows"][f"{kind}|{arm}|{enc}|{B}|{s}"][field]
                                  for s in SEEDS]) for p in P])

    print(f"\n=== D1  ENCODING vs STEP COUNT, n = {len(P)} targets, 4 seeds, "
          f"track_quality OFF ===")
    print("Basis: BUILT CHAIN (ideal-geometry CA trace from the arm's torsions). "
          "RMSD is ORACLE, post-hoc.\n")

    print("STEP COUNT actually realised (mean iterations):")
    print(f"  {'kind':<6}{'arm':<9}{'B':>6}" + "".join(f"{e:>10}" for e in ENCS) + "   ratio")
    for kind in KINDS:
        for arm in ARMS:
            for B in BUDGETS:
                v = [cell(kind, arm, e, B, "iters").mean() for e in ENCS]
                print(f"  {kind:<6}{arm:<9}{B:>6}" + "".join(f"{x:>10.1f}" for x in v)
                      + f"   {v[1]/max(v[0],1e-9):.2f}")

    print("\nDISPLACEMENT actually realised, ||wrap(z_best - z0)|| in radians (TORSION space,")
    print("the same basis for both encodings -- s20 recorded update_norm in each arm's OWN space):")
    print(f"  {'kind':<6}{'arm':<9}{'B':>6}" + "".join(f"{e:>10}" for e in ENCS) + "   ratio")
    for kind in KINDS:
        for arm in ARMS:
            for B in BUDGETS:
                v = [cell(kind, arm, e, B, "disp_best").mean() for e in ENCS]
                print(f"  {kind:<6}{arm:<9}{B:>6}" + "".join(f"{x:>10.3f}" for x in v)
                      + f"   {v[1]/max(v[0],1e-9):.2f}")

    print("\nGAUGE RADIUS DRIFT ||u|| (embedded only; 1.000 at the start, and the objective's")
    print("gradient has EXACTLY zero radial component, so any drift is pure noise injection):")
    print(f"  {'kind':<6}{'arm':<9}{'B':>6}{'r_last':>10}{'r_max':>10}")
    for kind in KINDS:
        for arm in ARMS:
            for B in BUDGETS:
                print(f"  {kind:<6}{arm:<9}{B:>6}"
                      f"{cell(kind,arm,'emb',B,'r_last').mean():>10.3f}"
                      f"{cell(kind,arm,'emb',B,'r_last_max').mean():>10.3f}")

    print("\nRMSD by encoding and budget (mean over targets and 4 seeds):")
    print(f"  {'kind':<6}{'arm':<9}" + "".join(f"{e+'@'+str(B):>13}" for e in ENCS
                                               for B in BUDGETS))
    for kind in KINDS:
        for arm in ARMS:
            print(f"  {kind:<6}{arm:<9}" + "".join(f"{cell(kind,arm,e,B).mean():>13.3f}"
                                                   for e in ENCS for B in BUDGETS))

    print("\n--- A: SPRINT 20's COMPARISON, matched BUDGET (emb@512 - theta@512) ---")
    A = {}
    for kind in KINDS:
        for arm in ARMS:
            dd = cell(kind, arm, "emb", 512) - cell(kind, arm, "theta", 512)
            mu, lo, hi = _boot(dd, rng)
            A[(kind, arm)] = dd
            print(f"  {kind:<6}{arm:<9}{mu:+7.3f} [{lo:+.3f},{hi:+.3f}]  "
                  f"W/L {(dd<0).sum()}/{(dd>0).sum()}  med {np.median(dd):+.3f}  "
                  f"folds {len({int(f) for f in np.unique(folds)})}")

    print("\n--- B: THE OPERATIVE CONTROL, matched STEP COUNT, both directions ---")
    print("   B1  theta handicapped down to the embedding's step count")
    print("   B2  the embedding given the budget that matches theta's step count")
    out = {"n": len(P), "targets": P, "cfg": CFG, "A": {}, "B1": {}, "B2": {}}
    for kind in KINDS:
        for arm in ARMS:
            it_e = cell(kind, arm, "emb", 512, "iters").mean()
            # pick the theta budget whose step count is closest to emb@512's
            cand = {B: abs(cell(kind, arm, "theta", B, "iters").mean() - it_e) for B in BUDGETS}
            Bt = min(cand, key=cand.get)
            it_t = cell(kind, arm, "theta", 512, "iters").mean()
            cand2 = {B: abs(cell(kind, arm, "emb", B, "iters").mean() - it_t) for B in BUDGETS}
            Be = min(cand2, key=cand2.get)
            d1 = cell(kind, arm, "emb", 512) - cell(kind, arm, "theta", Bt)
            d2 = cell(kind, arm, "emb", Be) - cell(kind, arm, "theta", 512)
            m1 = _boot(d1, rng); m2 = _boot(d2, rng)
            mA = _boot(A[(kind, arm)], rng)
            out["A"][f"{kind}|{arm}"] = {"mean": mA[0], "ci95": [mA[1], mA[2]],
                                         "W": int((A[(kind,arm)]<0).sum()),
                                         "L": int((A[(kind,arm)]>0).sum())}
            out["B1"][f"{kind}|{arm}"] = {"theta_budget": Bt, "emb_iters": it_e,
                                          "theta_iters": cell(kind,arm,"theta",Bt,"iters").mean(),
                                          "mean": m1[0], "ci95": [m1[1], m1[2]],
                                          "W": int((d1<0).sum()), "L": int((d1>0).sum())}
            out["B2"][f"{kind}|{arm}"] = {"emb_budget": Be, "theta_iters": it_t,
                                          "emb_iters": cell(kind,arm,"emb",Be,"iters").mean(),
                                          "mean": m2[0], "ci95": [m2[1], m2[2]],
                                          "W": int((d2<0).sum()), "L": int((d2>0).sum())}
            print(f"  {kind:<6}{arm:<9}"
                  f"B1 emb@512 vs theta@{Bt} ({it_e:.0f} vs "
                  f"{cell(kind,arm,'theta',Bt,'iters').mean():.0f} steps) "
                  f"{m1[0]:+7.3f} [{m1[1]:+.3f},{m1[2]:+.3f}] W/L {(d1<0).sum()}/{(d1>0).sum()}")
            print(f"  {'':<6}{'':<9}"
                  f"B2 emb@{Be} vs theta@512 ({cell(kind,arm,'emb',Be,'iters').mean():.0f} vs "
                  f"{it_t:.0f} steps) "
                  f"{m2[0]:+7.3f} [{m2[1]:+.3f},{m2[2]:+.3f}] W/L {(d2<0).sum()}/{(d2>0).sum()}")

    print("\n--- C: THE STEP-COUNT RESPONSE ITSELF, within each encoding "
          "(does optimising harder hurt?) ---")
    out["C"] = {}
    if len(BUDGETS) < 2:
        out["C"] = {"status": "NOT RUN -- single budget in this configuration"}
        print("  (single budget in this configuration -- the step-count response is NOT RUN,")
        print("   and it is recorded as NOT RUN rather than omitted)")
    for kind in (KINDS if len(BUDGETS) >= 2 else ()):
        for arm in ARMS:
            for e in ENCS:
                hiB, loB = max(BUDGETS), min(BUDGETS)
                dd = cell(kind, arm, e, hiB) - cell(kind, arm, e, loB)
                mu, lo, hi = _boot(dd, rng)
                out["C"][f"{kind}|{arm}|{e}"] = {"mean": mu, "ci95": [lo, hi], "hiB": hiB,
                                                 "loB": loB,
                                                 "W": int((dd<0).sum()), "L": int((dd>0).sum())}
                print(f"  {kind:<6}{arm:<9}{e:<7}B={hiB} minus B={loB}: {mu:+7.3f} "
                      f"[{lo:+.3f},{hi:+.3f}]  W/L {(dd<0).sum()}/{(dd>0).sum()}")

    print("\n--- D: TARGET-LEVEL POOLED, the brief's mandated unit ---")
    for tag, get in (("A matched budget", lambda k, a: A[(k, a)]),
                     ("B1 matched steps", lambda k, a: cell(k, a, "emb", 512)
                      - cell(k, a, "theta", out["B1"][f"{k}|{a}"]["theta_budget"])),
                     ("B2 matched steps", lambda k, a: cell(k, a, "emb",
                      out["B2"][f"{k}|{a}"]["emb_budget"]) - cell(k, a, "theta", 512))):
        M = np.stack([get(k, a) for k in KINDS for a in ARMS])
        per_t = M.mean(0)
        mu, lo, hi = _boot(per_t, rng)
        out.setdefault("pooled", {})[tag] = {"mean": mu, "ci95": [lo, hi],
                                             "W": int((per_t<0).sum()),
                                             "L": int((per_t>0).sum()),
                                             "sign_p": _signp(len(per_t), int((per_t>0).sum()))}
        print(f"  {tag:<20}{mu:+7.4f} [{lo:+.4f},{hi:+.4f}]  W/L "
              f"{(per_t<0).sum()}/{(per_t>0).sum()}  sign p = "
              f"{_signp(len(per_t), int((per_t>0).sum())):.4f}")

    print("\nMDE at 80% power = 0.084 A.  A residual inside it is NOT MEASURED, never 'matched'.")
    json.dump(out, open(os.path.join(RESULTS, f"{TAG}_report.json"), "w"), indent=1)
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if len(sys.argv) > 5:                       # e.g. ... run LEG,DIST spsa,adam_fd 512 60
        #: POWER.  Sec.4.5 of the findings shows an n=10 panel's own 80%-power MDE on this
        #: comparison is 0.24-0.80 A -- it cannot resolve the encoding question at all, and
        #: neither could Sprint 20's.  On the objectives that are cheap enough, use the tuning
        #: instrument's first N targets in pinned pdb-sorted order instead of the 10-target
        #: panel.  This is no longer paired to Sprint 20's panel; it is an independent and
        #: better-powered replication, and it is labelled as one.
        from s12 import instrument as _I
        TARGETS = [t["pdb"] for t in _I.targets()][:int(sys.argv[5])]
        CFG["targets"] = TARGETS
        GATE_TARGETS = TARGETS[:3]
    if len(sys.argv) > 3:                       # e.g. ... run AMBc spsa 512
        ARMS = tuple(sys.argv[3].split(","))
        CFG["arms"] = list(ARMS)
    if len(sys.argv) > 4:
        BUDGETS = tuple(int(x) for x in sys.argv[4].split(","))
        CFG["budgets"] = list(BUDGETS)
    if len(sys.argv) > 2:                       # e.g.  python -m s21.d_enc run LEG
        KINDS = tuple(sys.argv[2].split(","))
        CFG["kinds"] = list(KINDS)
        N_CELLS = (len(TARGETS) * len(KINDS) * len(ARMS) * len(ENCS)
                   * len(BUDGETS) * len(SEEDS))
        TAG = ("d_enc_" + "-".join(KINDS)
               + ("" if len(sys.argv) <= 5 else f"_n{len(TARGETS)}")
               + ("" if len(sys.argv) <= 3 else "_" + "-".join(ARMS))
               + ("" if len(sys.argv) <= 4 else "_B" + "-".join(str(b) for b in BUDGETS)))
    if mode == "gate":
        gate()
    elif mode == "report":
        report()
    else:
        run()
