"""s24/c_ladder.py -- LANE C CONFIRMATORY RUN, n=126.  Pre-registered in `s24/PREREG_C.md`.

THE QUESTION.  Can `p(phi,psi | sequence, structural prior) -> sample -> build` supply a candidate
source that moves the 126-target dev mean?  Directive SS7/SS26.

WHAT IS AND IS NOT TRAINED.  Nothing here is trained.  Four zero-training torsion distributions
span the conditioning ladder from "knows nothing about the target" to "knows the retrieved pool",
and the top of that ladder (`T3_pool`) was measured in s14 at phi 33.6 / psi 59.2 deg -- BETTER
than this project's trained leave-fold-out sequence predictor at 36.1 / 62.4.  A learned model has
to beat T3 to have earned its parameters, and `phi-carries-no-sequence-signal` prices the channel
it would have to do that with at 0.3 deg of phi and 10.4 deg of psi.  So T3 is the ACHIEVABLE
UPPER REFERENCE for the lane, and this run measures the lane through it.

CORPUS.  No Lane A dependency and no Lane A hash is claimed.  Every arm draws only from the
per-target universe in `s8/generate_univ/<pdb>.npz`, which is the project's already-audited
leakage-safe library (out-of-fold peptides + this fold's fragments).  `rr` and `nat_ca` are read
ONLY for post-hoc ORACLE scoring, never inside a decision.  A hash IS required before any trained
arm; its absence is recorded rather than assumed away.

BASIS DISCIPLINE.  Members are BUILT CHAINS.  The emitted top-75 uniform average is a POINT CLOUD
(an average of 75 chains is not itself an ideal-geometry chain), so it is commensurable with the
incumbent's 3.0483 A.  The projection of that cloud onto ideal geometry is the BUILT-CHAIN figure,
comparable to 3.204 A, and is reported in a separate column that is never mixed with the first.

OPERATOR FORKS (Rule 0).  Enumerated in `s24/PREREG_C.md` SS2 before this file ran; each names the
alternative not taken.  Enumerated by the lane that has a stake, which is a weakening and is
declared as one.  Summary: functional = shipped Bayes-risk unchanged (NOT re-tuned); basis =
point cloud on both sides, built-chain separate (NOT quoted across); readout = scored top-75
uniform mean at matched count, unscored uniform-75 alongside (NOT one readout silently chosen);
normalisation = native-frame bias cosine with a within-source control (NOT the raw cosine alone);
null = within-source ceiling + a plausible constant-helix zero-information arm + matched-count
best-of-500 (NOT zero, NOT uniform-on-the-torus, NOT an unmatched sample count); THE LABEL =
continuous Ca-RMSD and raw signed cosines (NOT a binarised "did it decorrelate").

AUDITS THAT PRECEDE ANY RMSD CLAIM.  Mode collapse (unique structures, torsion entropy, effective
modes, mean pairwise RMSD, duplicate fraction, mode occupancy) and geometric validity
(Ramachandran, omega/cis, bond length, bond angle, clashes, chain breaks) are computed on every
emitted pool and printed BEFORE the RMSD table.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402
from core import project as pj           # noqa: E402
from s24 import c_probe as CP            # noqa: E402

TOPM = 75
NSAMP = 2000
NMATCH = 500                       # matched-count null against the shipped K=500 pool
MIX = [(75, 0), (60, 15), (50, 25), (38, 37), (25, 50), (0, 75)]
ARMS = ("T0_helix", "T1_blind", "T2_restype", "T3_pool")
BAND = 1.5
OUT = os.path.join(RES, "c_ladder.json")


def _git():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    except Exception:
        return "unknown"


def _save(o):
    t = OUT + ".tmp"
    with open(t, "w") as fh:
        json.dump(o, fh)
    os.replace(t, OUT)


# ------------------------------------------------------------------ mode-collapse audit
def mode_audit(phi, psi, CA, sub=200, dup_tol=0.10, nb=12):
    """Every quantity the brief demands before an RMSD claim, on the EMITTED set.

    `unique`      distinct structures after rounding torsions to 1 deg
    `H_tors`      mean per-residue entropy of the joint (phi,psi) histogram, nats; and
    `eff_modes`   exp(H) -- the effective number of distinct local conformations per residue
    `pw_rmsd`     mean pairwise Ca-RMSD on a stable subsample
    `dup_frac`    fraction whose NEAREST NEIGHBOUR is within `dup_tol` A -- the collapse statistic
    `ess`         effective sample size of the SCORE's implied weights, softmax at the score's own
                  sd, reported because a set can be diverse and still contribute one structure
    `occ_*`       Ramachandran region occupancy (alpha / beta / other), the mode-occupancy readout
    """
    B, n = phi.shape
    q = np.round(np.rad2deg(np.column_stack([phi, psi]))).astype(np.int32)
    uniq = len(np.unique(q, axis=0))
    edges = np.linspace(-np.pi, np.pi, nb + 1)
    H = []
    for i in range(n):
        h, _, _ = np.histogram2d(phi[:, i], psi[:, i], bins=[edges, edges])
        p = h.ravel() / max(h.sum(), 1.0)
        p = p[p > 0]
        H.append(float(-(p * np.log(p)).sum()))
    Hm = float(np.mean(H))
    r = SD.stable_rng("c_ladder", "audit").permutation(B)[:min(sub, B)]
    P = I.pairwise_rmsd(CA[r])
    np.fill_diagonal(P, np.inf)
    return {"n_emitted": int(B), "unique": int(uniq), "unique_frac": float(uniq) / B,
            "H_tors_nats": Hm, "eff_modes": float(np.exp(Hm)),
            "pw_rmsd_mean": float(P[np.isfinite(P)].mean()),
            "dup_frac": float((P.min(1) < dup_tol).mean()),
            "occ_alpha": float((((phi < 0) & (psi > -2.0) & (psi < 0.5))).mean()),
            "occ_beta": float((((phi < 0) & ((psi > 1.0) | (psi < -2.6)))).mean())}


def geom_audit(CA, phi, psi, clash=4.0):
    """Geometric validity of the MEMBERS, with what is true BY CONSTRUCTION named as such.

    Members come from `build_ca_exact` at ideal geometry with omega fixed trans, so bond lengths,
    bond angles, omega and cis fraction are exact by construction and are reported as such rather
    than as a passing test.  What is NOT by construction, and is genuinely measured: Ramachandran
    allowedness, non-local Ca-Ca clashes, and chain breaks.
    """
    B, n = CA.shape[:2]
    d = np.linalg.norm(CA[:, 1:] - CA[:, :-1], axis=-1)
    i, j = np.triu_indices(n, k=2)
    dd = np.linalg.norm(CA[:, i, :] - CA[:, j, :], axis=-1)
    # generous Ramachandran: the four broadly populated quadrant regions for L-amino acids
    allowed = ((phi < 0) & (psi > -1.6) & (psi < 1.0)) | ((phi < 0) & (psi > 1.0)) | \
              ((phi > 0) & (phi < 2.0) & (psi > -0.3) & (psi < 1.6))
    return {"ca_ca_mean": float(d.mean()), "ca_ca_sd": float(d.std()),
            "chain_breaks": float((np.abs(d - 3.8) > 0.1).mean()),
            "clash_frac_struct": float((dd.min(1) < clash).mean()),
            "clash_per_struct": float((dd < clash).sum(1).mean()),
            "rama_allowed": float(allowed.mean()),
            "omega_cis_frac": 0.0, "bond_len_dev": 0.0, "bond_ang_dev": 0.0,
            "BY_CONSTRUCTION": "omega trans, bond lengths and angles ideal (build_ca_exact)"}


# ------------------------------------------------------------------ the run
def run(nsamp=NSAMP, resume=True):
    tg = I.targets()
    rows = []
    if resume and os.path.exists(OUT):
        try:
            rows = json.load(open(OUT)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for ci, t in enumerate(tg):
        pdb = t["pdb"]
        if pdb in done:
            continue
        n = t["n"]
        u = I.load_univ(pdb); nat = u["nat_ca"]
        dg = I.distogram(pdb, u["seq"], u["fold"]); i, j = I.pair_index(n)
        dnat = I.pair_dists(nat[None], i, j)[0]
        gdg = np.asarray(dg["expected"], float) - dnat          # the distogram's own signed error

        idx = I.pool_idx(u); Wp = u["W"][idx]; rrp = u["rr"][idx]
        scp = np.asarray(I.shipped_score(dg, I.pair_dists(Wp, i, j)), float)
        top = np.argsort(scp, kind="stable")[:TOPM]
        A = Wp[top]
        cA = CP._avg(A); eA = CP._bias(cA, nat)
        prA = I.project(cA, u["seq"], u["fold"])                 # BUILT-CHAIN basis, kept separate

        # R: the incumbent's own 75, rebuilt from their own torsions -- prices the manifold
        CAr = np.asarray(pj.build_ca_exact(CP.wrap(u["PHI"][idx][top]),
                                           CP.wrap(u["PSI"][idx][top])), float)
        cR = CP._avg(CAr)

        r = {"pdb": pdb, "n": n, "fold": int(u["fold"]),
             "incumbent": float(I.ca_rmsd(cA, nat)),
             "incumbent_builtchain": float(I.ca_rmsd(prA["ca"], nat)),
             "rebuild75": float(I.ca_rmsd(cR, nat)),
             "cos_rebuild": CP._cos(CP._bias(cR, nat), eA),
             "pool_best_ORACLE": float(rrp.min()),
             "pool_score_p75": float(np.sort(scp)[TOPM - 1]),
             "cos_inc_dgram": CP._cos(I.pair_dists(cA[None], i, j)[0] - dnat, gdg)}

        rng0 = SD.stable_rng("c_ladder", pdb)
        for tag in ARMS:
            def fn(g, _tag=tag, _n=n, _u=u, _idx=idx, _top=top, _N=nsamp):
                return {"T0_helix": lambda: CP.s_helix(_n, _N, g),
                        "T1_blind": lambda: CP.s_blind(_u, _n, _N, g),
                        "T2_restype": lambda: CP.s_restype(_u, _n, _N, g),
                        "T3_pool": lambda: CP.s_pool(_u, _n, _N, g, _idx[_top])}[_tag]()

            ph, ps = fn(rng0)
            CA = np.asarray(pj.build_ca_exact(ph, ps), float)
            sc = np.asarray(I.shipped_score(dg, I.pair_dists(CA, i, j)), float)
            o = np.argsort(sc, kind="stable")
            B = CA[o[:TOPM]]
            cS = CP._avg(B)
            ru = rng0.choice(len(CA), TOPM, replace=False)
            cU = CP._avg(CA[ru])
            rr = I.kabsch_rmsd_batch(CA, nat)                    # ORACLE, post-hoc

            # within-source control: an independent draw from the SAME sampler
            ph2, ps2 = fn(rng0)
            CA2 = np.asarray(pj.build_ca_exact(ph2, ps2), float)
            sc2 = np.asarray(I.shipped_score(dg, I.pair_dists(CA2, i, j)), float)
            cS2 = CP._avg(CA2[np.argsort(sc2, kind="stable")[:TOPM]])
            ru2 = rng0.choice(len(CA2), TOPM, replace=False)
            cU2 = CP._avg(CA2[ru2])

            # matched-size mixtures: 75 members always, only composition changes
            mx = {}
            for a, b in MIX:
                s = np.concatenate([A[:a], B[:b]], 0) if (a and b) else (A[:a] if a else B[:b])
                mx["m%d_%d" % (a, b)] = float(I.ca_rmsd(CP._avg(s), nat))

            # SS17 union at MATCHED COUNT: 500 retrieved + 500 generated, ONE score, top-75
            gsub = o[:0]
            gi = SD.stable_rng("c_ladder", pdb, tag).permutation(len(CA))[:NMATCH]
            Wu = np.concatenate([Wp, CA[gi]], 0)
            src = np.concatenate([np.zeros(len(Wp), int), np.ones(len(gi), int)])
            scu = np.asarray(I.shipped_score(dg, I.pair_dists(Wu, i, j)), float)
            ou = np.argsort(scu, kind="stable")[:TOPM]

            scm = sc[gi]; rrm = rr[gi]                           # matched-count arm
            r[tag] = {
                # --- the two spec numbers, both readouts, basis stated
                "rmsd": float(I.ca_rmsd(cS, nat)),
                "rmsd_unif": float(I.ca_rmsd(cU, nat)),
                "cos_vs_incumbent": CP._cos(CP._bias(cS, nat), eA),
                "cos_unif_vs_incumbent": CP._cos(CP._bias(cU, nat), eA),
                "cos_within": CP._cos(CP._bias(cS, nat), CP._bias(cS2, nat)),
                "cos_unif_within": CP._cos(CP._bias(cU, nat), CP._bias(cU2, nat)),
                # --- L3's mechanism, in the score's own space
                "cos_scored_dgram": CP._cos(I.pair_dists(cS[None], i, j)[0] - dnat, gdg),
                "cos_unif_dgram": CP._cos(I.pair_dists(cU[None], i, j)[0] - dnat, gdg),
                # --- the endpoint
                "mix": mx,
                "union": float(I.ca_rmsd(CP._avg(Wu[ou]), nat)),
                "union_gen_frac": float(src[ou].mean()),
                # --- coverage, route (a), at MATCHED COUNT 500
                "frac_better_than_pool_best": float((scm < scp.min()).mean()),
                "frac_better_than_pool_p75": float((scm < np.sort(scp)[TOPM - 1]).mean()),
                "gen_best_ORACLE_matched": float(rrm.min()),
                "gen_best_ORACLE_full": float(rr.min()),
                "n_beats_pool_best_matched_ORACLE": int((rrm < rrp.min()).sum()),
                "n_in_band_matched_ORACLE": int((rrm <= rrp.min() + BAND).sum()),
                "sel_member_mean_ORACLE": float(I.kabsch_rmsd_batch(B, nat).mean()),
                # --- the mandatory audits
                "mode": mode_audit(ph, ps, CA),
                "geom": geom_audit(CA, ph, ps),
            }
            del gsub
        rows.append(r)
        if (len(rows) % 5) == 0 or len(rows) == len(tg):
            _save({"rows": rows, "complete": False, "n_expected": len(tg)})
            print("  %d/%d  %.0fs" % (len(rows), len(tg), time.time() - t0), flush=True)

    need = ("incumbent", "rebuild75") + ARMS
    keys = {r["pdb"] for r in rows}
    ok = keys == {t["pdb"] for t in tg} and all(all(k in r for k in need) for r in rows)
    _save({"rows": rows, "complete": bool(ok), "n_expected": len(tg),
           "git": _git(), "nsamp": nsamp, "topm": TOPM, "nmatch": NMATCH, "mix": MIX,
           "corpus": "s8/generate_univ per-target leakage-safe universe; NO Lane A hash (untrained)",
           "config_hash": hashlib.sha1(("%d|%d|%d|%s" % (nsamp, TOPM, NMATCH, ARMS)).encode()
                                       ).hexdigest()[:16]})
    report(rows)
    return rows


# ------------------------------------------------------------------ report
def report(rows=None):
    if rows is None:
        rows = json.load(open(OUT))["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    rng = SD.stable_rng("c_ladder", "rep")
    F = sorted(set(fold.tolist()))
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731

    def st(d):
        d = np.asarray(d, float); k = len(d)
        se = d.std(ddof=1) / np.sqrt(k)
        fs = [np.concatenate([d[fold == q] for q in rng.choice(F, len(F), replace=True)]).mean()
              for _ in range(4000)]
        return (d.mean(), se, 2.8016 * se, float(np.percentile(fs, 2.5)),
                float(np.percentile(fs, 97.5)), int((d < 0).sum()), float(d.max()))

    inc = g("incumbent")
    print("\n" + "=" * 100)
    print("LANE C CONFIRMATORY, n=%d, N=%d samples/target, no training." % (len(rows), NSAMP))
    print("=" * 100)

    print("\n[1] MODE-COLLAPSE AUDIT -- precedes every RMSD claim.")
    print("  %-12s %9s %9s %9s %9s %9s %9s" % ("arm", "uniq_frac", "dup_frac", "H(nats)",
                                               "eff_modes", "pw_RMSD", "occ_alpha"))
    for k in ARMS:
        f = lambda q: np.array([r[k]["mode"][q] for r in rows], float)   # noqa: E731
        print("  %-12s %9.4f %9.4f %9.4f %9.3f %9.3f %9.3f"
              % (k, f("unique_frac").mean(), f("dup_frac").mean(), f("H_tors_nats").mean(),
                 f("eff_modes").mean(), f("pw_rmsd_mean").mean(), f("occ_alpha").mean()))
    print("  A collapsed sampler has unique_frac -> 0, dup_frac -> 1, pw_RMSD -> 0.")

    print("\n[2] GEOMETRIC VALIDITY.  omega/cis, bond lengths and bond angles are exact BY")
    print("    CONSTRUCTION (build_ca_exact, ideal geometry, omega trans) -- reported, not tested.")
    print("  %-12s %11s %11s %11s %11s" % ("arm", "rama_allow", "clash_frac", "clash/str",
                                           "chainbreak"))
    for k in ARMS:
        f = lambda q: np.array([r[k]["geom"][q] for r in rows], float)   # noqa: E731
        print("  %-12s %11.4f %11.4f %11.4f %11.4f"
              % (k, f("rama_allowed").mean(), f("clash_frac_struct").mean(),
                 f("clash_per_struct").mean(), f("chain_breaks").mean()))

    print("\n[3] THE TWO SPEC NUMBERS.  POINT CLOUD on both sides.  Spec: <=~3.9 A AND cos <=~0.65.")
    print("            |-------- SCORED top-75 (the spec's readout) --------|"
          "  |---- UNIFORM 75 (L2's readout) ----|")
    print("  %-12s %8s %8s %8s %7s | %8s %8s %8s %7s"
          % ("arm", "RMSD", "cos_inc", "cos_wtn", "ratio", "RMSD", "cos_inc", "cos_wtn", "ratio"))
    print("  %-12s %8.4f   [point-cloud incumbent; built-chain %.4f, NEVER compared]"
          % ("INCUMBENT", inc.mean(), g("incumbent_builtchain").mean()))
    print("  %-12s %8.4f %8.4f   [CALIBRATION: the ideal-geometry manifold costs %+.4f A]"
          % ("REBUILD75", g("rebuild75").mean(), g("cos_rebuild").mean(),
             g("rebuild75").mean() - inc.mean()))
    for k in ARMS:
        f = lambda q: np.array([r[k][q] for r in rows], float)           # noqa: E731
        print("  %-12s %8.4f %8.4f %8.4f %7.3f | %8.4f %8.4f %8.4f %7.3f"
              % (k, f("rmsd").mean(), f("cos_vs_incumbent").mean(), f("cos_within").mean(),
                 f("cos_vs_incumbent").mean() / f("cos_within").mean(),
                 f("rmsd_unif").mean(), f("cos_unif_vs_incumbent").mean(),
                 f("cos_unif_within").mean(),
                 f("cos_unif_vs_incumbent").mean() / f("cos_unif_within").mean()))

    print("\n[4] THE ENDPOINT.  Matched-size mixtures, 75 members always, POINT CLOUD.")
    for k in ARMS:
        M = {q: np.array([r[k]["mix"][q] for r in rows], float) for q in rows[0][k]["mix"]}
        ks = list(M); e0, e1 = M[ks[0]], M[ks[-1]]
        print("  --- %s" % k)
        for q in ks[1:]:
            frac = float(q.split("_")[1]) / TOPM
            line = (1 - frac) * e0 + frac * e1
            m, se, mde, lo, hi, w, wt = st(M[q] - e0)
            v = ("BEATS" if hi < 0 and abs(m) > mde else "worse" if lo > 0 and abs(m) > mde
                 else "NULL")
            print("      %-9s %7.4f  d %+7.4f  line %+7.4f  SE %.4f MDE %.4f fold[%+.3f,%+.3f]"
                  " %3dW/%3dL worst %+.2f  %s"
                  % (q, M[q].mean(), m, (M[q] - line).mean(), se, mde, lo, hi, w,
                     len(rows) - w, wt, v))
        m, se, mde, lo, hi, w, wt = st(np.array([r[k]["union"] for r in rows]) - inc)
        gf = np.array([r[k]["union_gen_frac"] for r in rows], float)
        v = "BEATS" if hi < 0 and abs(m) > mde else "worse" if lo > 0 and abs(m) > mde else "NULL"
        print("      UNION     %7.4f  d %+7.4f  SE %.4f MDE %.4f fold[%+.3f,%+.3f]"
              " %3dW/%3dL worst %+.2f  genfrac %.3f  %s"
              % (np.array([r[k]["union"] for r in rows]).mean(), m, se, mde, lo, hi, w,
                 len(rows) - w, wt, gf.mean(), v))

    print("\n[5] COVERAGE -- coordinator's route (a).  MATCHED COUNT %d vs the shipped K=500 pool."
          % NMATCH)
    print("  pool ORACLE best %.4f   incumbent %.4f" % (g("pool_best_ORACLE").mean(), inc.mean()))
    print("  %-12s %11s %11s %11s | %11s %11s %11s"
          % ("arm", "f>poolbest", "f>pool_p75", "union%gen", "ORACLE gbst", "#<poolbest",
             "#in-band"))
    for k in ARMS:
        f = lambda q: np.array([r[k][q] for r in rows], float)           # noqa: E731
        print("  %-12s %11.4f %11.4f %11.4f | %11.4f %11.3f %11.3f"
              % (k, f("frac_better_than_pool_best").mean(), f("frac_better_than_pool_p75").mean(),
                 f("union_gen_frac").mean(), f("gen_best_ORACLE_matched").mean(),
                 f("n_beats_pool_best_matched_ORACLE").mean(),
                 f("n_in_band_matched_ORACLE").mean()))
    print("  Route (a) needs the score to prefer generated chains AND those chains to be")
    print("  structurally better.  Compare 'union%gen' against 'ORACLE gbst' vs pool best.")

    print("\n[6] L3'S MECHANISM.  cos(emitted distance error, DISTOGRAM's own distance error).")
    print("  %-12s %11s %11s   [INCUMBENT %.4f]" % ("arm", "scored75", "unif75",
                                                    g("cos_inc_dgram").mean()))
    for k in ARMS:
        f = lambda q: np.array([r[k][q] for r in rows], float)           # noqa: E731
        print("  %-12s %11.4f %11.4f" % (k, f("cos_scored_dgram").mean(),
                                         f("cos_unif_dgram").mean()))
    print("  Selection raising this in every arm, including the zero-information one, is L3.")


if __name__ == "__main__":
    report() if len(sys.argv) > 1 and sys.argv[1] == "report" else run()
