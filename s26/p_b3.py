"""s26/p_b3.py -- PREREG_B3: where does the pipeline gain over sequence-only, and is that set
characterisable native-free?

Arms (all BUILT CHAINS, persisted; nothing is re-emitted):
    arm    production rmsd_arm  (bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json)      3.2148
    tors   sequence-only torsion predictor, direct build (s13/cache/tors_rows.npz a_pepPos) 3.7705
    helix  constant alpha-helix, zero information (s14/results/ladder.json L0)           4.0648
Label: d = arm - tors (negative = the pipeline is better; 89/126), sign(d), and d_helix.

NATIVE-FREE FEATURES (built by `features`, cached to s26/results/p_b3_features.json; allowed
before sign-off because nothing in them reads a native): length; composition (20 fractions,
hydrophobic / charged / Gly+Pro fractions, net charge); the retrieval pool's secondary-structure
class (H/E/C fractions of the shipped top-75 members' own torsions through I.ss_of, and the
share of members whose majority class is H, E, C); the shipped distogram's per-pair bin entropy
(mean, max, sd, multimodal fraction, mean sd, expected mean distance / n); ESM-2 650M contact-map
statistics (mean, long-range mass, max, entropy); retrieval-score entropy and normalised gaps.

TEST (PREREG_B3): nested leave-fold-out ridge classifier of sign(d) against a label-permutation
null (balanced accuracy vs the null's 95th percentile), and nested ridge regression of d (held-out
R^2; per-target MSE reduction with the fold-clustered CI through ST.compare).  `run` reads the
persisted RMSDs and is gated on "PHASE 0 SIGNED OFF".

    python s26/p_b3.py features      (native-free; allowed now)
    python s26/p_b3.py selftest      (synthetic labels on the REAL features; no RMSD read)
    python s26/p_b3.py run           (after sign-off)
"""
from __future__ import annotations

import glob
import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from s26 import p_stats as PS                        # noqa: E402
from s26 import p_ladder as L                        # noqa: E402

RES = os.path.join(HERE, "results")
FEAT = os.path.join(RES, "p_b3_features.json")
OUT = os.path.join(RES, "p_b3.json")
AA = "ARNDCQEGHILKMFPSTWYV"
HYD, CHG, GP = set("AILMFVWY"), set("DEKR"), set("GP")
N_PERM = 300


def _entropy(p, axis=-1):
    p = np.asarray(p, float)
    return -(p * np.log(np.clip(p, 1e-12, 1.0))).sum(axis)


def _peaks(p, thr=0.02):
    p = np.asarray(p, float)
    left = np.r_[-np.inf, p[:-1]]; right = np.r_[p[1:], -np.inf]
    return int(((p > left) & (p >= right) & (p > thr)).sum())


def features_one(t, u, rec, dg, con):
    """One target's native-free feature dict.  `u` is used for W (pool), PHI/PSI, sim; never rr/nat_ca."""
    seq, n = t["seq"], int(t["n"])
    f = {"n": float(n)}
    for a in AA:
        f["aa_" + a] = seq.count(a) / n
    f["hyd_frac"] = sum(c in HYD for c in seq) / n
    f["chg_frac"] = sum(c in CHG for c in seq) / n
    f["gp_frac"] = sum(c in GP for c in seq) / n
    f["net_charge"] = (seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")) / n
    idx = I.pool_idx(u); sub = np.asarray(rec["sub"], int); mem = idx[sub]
    phi = np.asarray(u["PHI"][mem], float); psi = np.asarray(u["PSI"][mem], float)
    H = E = C = 0.0; majH = majE = majC = 0.0
    for k in range(len(mem)):
        ss = I.ss_of(phi[k], psi[k])
        h, e = ss.count("H") / n, ss.count("E") / n
        H += h; E += e; C += 1 - h - e
        mx = max(h, e, 1 - h - e)
        majH += mx == h; majE += mx == e and mx != h; majC += (mx == 1 - h - e) and mx != h and mx != e
    m = float(len(mem))
    f.update(ss_H=H / m, ss_E=E / m, ss_C=C / m, maj_H=majH / m, maj_E=majE / m, maj_C=majC / m)
    prob = np.asarray(dg["prob"], float)
    ent = _entropy(prob, 1)
    f.update(dg_ent_mean=float(ent.mean()), dg_ent_max=float(ent.max()), dg_ent_sd=float(ent.std()),
             dg_multimodal=float(np.mean([_peaks(p) >= 2 for p in prob])),
             dg_sd_mean=float(np.asarray(dg["sd"], float).mean()),
             dg_exp_over_n=float(np.asarray(dg["expected"], float).mean() / n))
    con = np.asarray(con, float); ii, jj = np.triu_indices(n, k=3)
    c3 = con[ii, jj]; lr = con[np.triu_indices(n, k=6)] if n > 6 else np.zeros(1)
    pc = c3 / max(c3.sum(), 1e-9)
    f.update(con_mean=float(c3.mean()), con_lr=float(lr.mean()), con_max=float(c3.max()),
             con_ent=float(_entropy(pc) / np.log(max(len(pc), 2))))
    sim = np.asarray(u["sim"][idx], float)
    z = sim - sim.max(); w = np.exp(z); w /= w.sum()
    f.update(sim_ent=float(_entropy(w) / np.log(len(w))), sim_gap01=float((sim[0] - sim[1]) / (sim.std() + 1e-9)),
             sim_gap75=float((sim[0] - sim[74]) / (sim.std() + 1e-9)), sim_mean_top75=float(sim[:75].mean()))
    return f


def build_features(verbose=True):
    if os.path.exists(FEAT):
        return json.load(open(FEAT))
    L.install_guard()
    con_t = L.contacts()
    rows = []
    for c, t in enumerate(I.targets()):
        u = I.load_univ(t["pdb"]); rec = I.shipped_record(t["pdb"]); dg = I.distogram(t["pdb"])
        f = features_one(t, u, rec, dg, con_t[t["seq"]])
        rows.append({"pdb": t["pdb"], "fold": int(t["fold"]), **f})
        del u
        if verbose and (c + 1) % 25 == 0:
            print("  features %d/126" % (c + 1), flush=True)
    names = [k for k in rows[0] if k not in ("pdb", "fold")]
    out = {"rows": rows, "names": names}
    ST.save_atomic(FEAT, out, complete_keys=names, rows=rows, n_expected=126, module_file=__file__)
    return out


def labels():
    """The persisted per-target arms.  Reads recorded RMSDs: call only after sign-off (or in selftest
    with synthetic labels, which never calls this)."""
    recs = {}
    for fpath in glob.glob(os.path.join(ROOT, "bench_results", "cache", I.PROD_KEY, "*.json")):
        r = json.load(open(fpath)); recs[r["pdb"]] = r
    z = np.load(os.path.join(ROOT, "s13", "cache", "tors_rows.npz"), allow_pickle=True)
    tors = {r["pdb"]: r for r in json.loads(str(z["a_pepPos"]))}
    hel = json.load(open(os.path.join(ROOT, "s14", "results", "ladder.json")))["per_target"]["L0_constant_helix"]
    return recs, tors, hel


def analyse(X, names, folds, d, y, pdbs, tag, rng, n_perm=N_PERM):
    """Classifier + regression with their nulls.  Pure function of arrays; used by run and selftest."""
    pc, a_c = PS.nested_predict(X, y, folds, "class")
    acc = PS.balanced_accuracy(y, pc)
    null = PS.permutation_null(X, y, folds, "class", n_perm, rng)
    pr, a_r = PS.nested_predict(X, d, folds, "reg")
    base = np.zeros_like(d)
    for f in sorted(set(folds.tolist())):
        base[folds == f] = d[folds != f].mean()
    se_model = (d - pr) ** 2; se_base = (d - base) ** 2
    r2 = 1 - se_model.sum() / se_base.sum()
    cmp = ST.compare(se_model, se_base, folds, names=pdbs, label="%s: held-out squared error, ridge - training-fold-mean" % tag)
    return {"tag": tag, "balanced_acc": acc, "null_p95": float(np.percentile(null, 95)),
            "null_mean": float(null.mean()), "p_perm": float((null >= acc).mean()),
            "alphas_class": a_c, "heldout_r2": float(r2), "alphas_reg": a_r,
            "mse_reduction": cmp, "n_perm": int(n_perm)}


def run(n_perm=N_PERM):
    if not L._signed_off():
        raise SystemExit("PHASE GATE: p_b3 run reads persisted RMSDs; refused before sign-off.")
    fz = build_features(); rows = fz["rows"]; names = fz["names"]
    recs, tors, hel = labels()
    pdbs = [r["pdb"] for r in rows]; folds = np.array([r["fold"] for r in rows])
    X = np.array([[r[k] for k in names] for r in rows], float)
    arm = np.array([recs[p]["rmsd_arm"] for p in pdbs]); tb = np.array([tors[p]["rmsd_build"] for p in pdbs])
    hx = np.array([hel[p] for p in pdbs])
    rng = np.random.default_rng(26)
    out = {"n": len(pdbs), "names": names, "means": {"arm": float(arm.mean()), "tors": float(tb.mean()), "helix": float(hx.mean())},
           "contrasts": {"tors_minus_arm": ST.compare(tb, arm, folds, names=pdbs, label="tors - arm (built chains)"),
                         "helix_minus_arm": ST.compare(hx, arm, folds, names=pdbs, label="helix - arm (built chains)")}}
    d = arm - tb; y = np.sign(d); y[y == 0] = 1
    out["vs_tors"] = analyse(X, names, folds, d, y, pdbs, "arm - tors", rng, n_perm)
    d2 = arm - hx; y2 = np.sign(d2); y2[y2 == 0] = 1
    out["vs_helix"] = analyse(X, names, folds, d2, y2, pdbs, "arm - helix", rng, n_perm)
    fail = set(I.FAIL18)
    out["oracle_stratum_FAIL18"] = {"n": int(sum(p in fail for p in pdbs)),
                                    "d_tors_FAIL18": float(np.mean([d[k] for k, p in enumerate(pdbs) if p in fail])),
                                    "d_tors_other": float(np.mean([d[k] for k, p in enumerate(pdbs) if p not in fail]))}
    ST.save_atomic(OUT, out, module_file=__file__)
    for k in ("vs_tors", "vs_helix"):
        r = out[k]
        print("%s: balanced acc %.3f vs null p95 %.3f (perm p %.3f); held-out R2 %+.3f" % (r["tag"], r["balanced_acc"], r["null_p95"], r["p_perm"], r["heldout_r2"]))
        print(ST.fmt(r["mse_reduction"]))
    return out


def selftest():
    """Synthetic labels on the REAL native-free features: a planted signal must be found, a random
    label must not.  No RMSD is read."""
    fz = build_features(verbose=False); rows = fz["rows"]; names = fz["names"]
    folds = np.array([r["fold"] for r in rows]); pdbs = [r["pdb"] for r in rows]
    X = np.array([[r[k] for k in names] for r in rows], float)
    rng = np.random.default_rng(1)
    Z = (X - X.mean(0)) / (X.std(0) + 1e-9)
    w = np.zeros(len(names)); w[names.index("ss_H")] = 1.0; w[names.index("n")] = 0.7
    d = Z @ w + 0.5 * rng.normal(size=len(pdbs)); y = np.sign(d)
    r = analyse(X, names, folds, d, y, pdbs, "planted", rng, n_perm=30)
    print("  planted: acc %.3f null p95 %.3f R2 %.3f" % (r["balanced_acc"], r["null_p95"], r["heldout_r2"]))
    assert r["balanced_acc"] > r["null_p95"] and r["heldout_r2"] > 0.5
    d = rng.normal(size=len(pdbs)); y = np.sign(d)
    r = analyse(X, names, folds, d, y, pdbs, "random", rng, n_perm=30)
    print("  random:  acc %.3f null p95 %.3f R2 %+.3f" % (r["balanced_acc"], r["null_p95"], r["heldout_r2"]))
    assert r["balanced_acc"] <= max(r["null_p95"] + 0.05, 0.62)
    print("  p_b3 selftest OK (%d features, %d targets)" % (len(names), len(pdbs)))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "features"
    {"features": lambda: build_features(), "selftest": selftest, "run": run}[cmd]()
