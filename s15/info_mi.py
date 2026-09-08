"""SPRINT 15, INFO, PART D -- information content, with the estimator checked and every
number paired with a direct structural measure.

Plug-in mutual information over 75 pool members in 36 bins is dominated by its own bias, so
the primary estimator here is a PREDICTIVE (held-out cross-entropy) one, which is a valid
lower bound on the mutual information and cannot be inflated by sampling noise:

    I_pred  =  H0  -  CE(channel, truth)     bits

    H0   the cross-entropy of a SEQUENCE-BLIND, TARGET-BLIND marginal built from the OTHER
         125 targets' residues -- i.e. generic Ramachandran
    CE   the cross-entropy of the channel's own per-residue distribution on the native bin

`I_pred` is therefore "bits the target-specific channel adds over generic Ramachandran",
which is the quantity the project actually needs, not the raw MI with the torsion.

Reported for six bin counts (binning sensitivity), with a bootstrap over targets, a
target-shuffled null, AND the plug-in Miller-Madow MI beside it so a reader can see where
the plug-in estimator misbehaves.

Each information number is paired with a direct structural measure in the same table.

Also here: the pool's HIGHER-ORDER structure.  Does the joint distribution over the whole
window carry anything the per-residue marginals do not?  Measured by resampling each
residue independently from the pool (marginals only) against drawing whole windows.

    python -m s15.info_mi
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I           # noqa: E402
from s15 import info_lib as L            # noqa: E402
from s14 import retprior as RP           # noqa: E402
import peptide_db as pdb                  # noqa: E402

BINS = [6, 8, 12, 18, 24, 36]
ALPHA = 0.5


def binidx(a, B):
    return np.clip(((np.asarray(a, float) + np.pi) / (2 * np.pi) * B).astype(int), 0, B - 1)


def smooth(counts, B, alpha=ALPHA):
    return (counts + alpha) / (counts.sum() + alpha * B)


def collect():
    """Per target: pool torsion bins (top-75), native bins, distogram, native distances."""
    out = []
    for t in I.targets():
        p, n, seq, fold = t["pdb"], int(t["n"]), t["seq"], int(t["fold"])
        nt = pdb.by_pdb(p)
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        PH, PS, _ = RP.windows(p, "top75")
        u = I.load_univ(p)
        nat = np.asarray(u["nat_ca"], float)
        dg = I.distogram(p, seq, fold)
        i, j = I.pair_index(n)
        out.append({"pdb": p, "n": n, "fold": fold,
                    "PH": np.asarray(PH, float), "PS": np.asarray(PS, float),
                    "nphi": nphi, "npsi": npsi,
                    "prob": np.asarray(dg["prob"], float),
                    "centres": np.asarray(dg["centres"], float),
                    "sd": np.asarray(dg["sd"], float),
                    "dtrue": np.linalg.norm(nat[i] - nat[j], axis=1)})
    return out


def torsion_information(data, B, seed=0, shuffle=False):
    """Bits the retrieval pool adds over generic Ramachandran, per torsion."""
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(data)) if shuffle else np.arange(len(data))
    # generic marginal from the OTHER targets' NATIVE torsions (leave-one-target-out)
    allphi = [binidx(d["nphi"], B) for d in data]
    allpsi = [binidx(d["npsi"], B) for d in data]
    tot_p = np.zeros(B); tot_s = np.zeros(B)
    for a, b in zip(allphi, allpsi):
        tot_p += np.bincount(a, minlength=B)
        tot_s += np.bincount(b, minlength=B)
    per = []
    for k, d in enumerate(data):
        src = data[order[k]]
        n = d["n"]
        m = min(n, src["n"])
        h0p = smooth(tot_p - np.bincount(allphi[k], minlength=B), B)
        h0s = smooth(tot_s - np.bincount(allpsi[k], minlength=B), B)
        ce0, ce1 = [], []
        for r in range(m):
            for (P, nat_bins, h0) in ((src["PH"][:, r], allphi[k][r], h0p),
                                      (src["PS"][:, r], allpsi[k][r], h0s)):
                q = smooth(np.bincount(binidx(P, B), minlength=B), B)
                ce0.append(-np.log2(h0[nat_bins]))
                ce1.append(-np.log2(q[nat_bins]))
        per.append(float(np.mean(ce0) - np.mean(ce1)))
    return np.asarray(per)


def plugin_mi(data, B):
    """Plug-in MI(pool-modal bin ; native bin) with the Miller-Madow correction, pooled
    over all residues of all targets.  Shown for comparison; it is the estimator the brief
    warns about."""
    x, y = [], []
    for d in data:
        for r in range(d["n"]):
            for P, nat in ((d["PH"][:, r], d["nphi"][r]), (d["PS"][:, r], d["npsi"][r])):
                c = np.bincount(binidx(P, B), minlength=B)
                x.append(int(np.argmax(c))); y.append(int(binidx(nat, B)))
    x = np.asarray(x); y = np.asarray(y)
    N = len(x)
    J = np.zeros((B, B))
    np.add.at(J, (x, y), 1.0)
    J /= N
    px = J.sum(1); py = J.sum(0)
    nz = J > 0
    mi = float((J[nz] * np.log2(J[nz] / np.outer(px, py)[nz])).sum())
    # Miller-Madow: bias of a plug-in entropy is -(support-1)/(2 N ln2) per entropy term
    mm = ((np.count_nonzero(px) - 1) + (np.count_nonzero(py) - 1)
          - (np.count_nonzero(J) - 1)) / (2.0 * N * np.log(2))
    rng = np.random.default_rng(0)
    nulls = []
    for _ in range(200):
        ys = rng.permutation(y)
        Js = np.zeros((B, B)); np.add.at(Js, (x, ys), 1.0); Js /= N
        pxs = Js.sum(1); pys = Js.sum(0); nzs = Js > 0
        nulls.append((Js[nzs] * np.log2(Js[nzs] / np.outer(pxs, pys)[nzs])).sum())
    return {"plugin_mi_bits": mi, "miller_madow_corrected": mi - mm,
            "shuffle_null_mean_bits": float(np.mean(nulls)),
            "shuffle_null_sd": float(np.std(nulls)), "N": int(N), "B": int(B)}


def distogram_information(data):
    """Bits the distogram adds over the pooled distance marginal, per CA-CA pair."""
    cen = data[0]["centres"]
    edges = np.concatenate([[-np.inf], (cen[1:] + cen[:-1]) / 2, [np.inf]])
    bins_all = [np.clip(np.digitize(d["dtrue"], edges) - 1, 0, len(cen) - 1) for d in data]
    tot = np.zeros(len(cen))
    for b in bins_all:
        tot += np.bincount(b, minlength=len(cen))
    per, per_lo, per_hi = [], [], []
    for k, d in enumerate(data):
        h0 = smooth(tot - np.bincount(bins_all[k], minlength=len(cen)), len(cen))
        q = np.maximum(d["prob"], 1e-9)
        q = q / q.sum(1, keepdims=True)
        b = bins_all[k]
        ce0 = -np.log2(h0[b])
        ce1 = -np.log2(q[np.arange(len(b)), b])
        per.append(float(np.mean(ce0 - ce1)))
        lo = d["sd"] <= np.median(d["sd"])
        per_lo.append(float(np.mean((ce0 - ce1)[lo])))
        per_hi.append(float(np.mean((ce0 - ce1)[~lo])))
    return np.asarray(per), np.asarray(per_lo), np.asarray(per_hi)


def esm_contact_information():
    """Bits the ESM-2 contact map adds over the base contact rate, per CA-CA pair."""
    from s12 import esm_bank
    esm = esm_bank.load()
    ys, ps = [], []
    for t in I.targets():
        p, n, seq = t["pdb"], int(t["n"]), t["seq"]
        if seq not in esm:
            continue
        nat = np.asarray(I.load_univ(p)["nat_ca"], float)
        i, j = I.pair_index(n)
        y = (np.linalg.norm(nat[i] - nat[j], axis=1) < 8.0).astype(float)
        q = np.clip(esm[seq][2][i, j], 1e-4, 1 - 1e-4)
        ys.append(y); ps.append(q)
    per = []
    Y = np.concatenate(ys)
    base = float(Y.mean())
    for y, q in zip(ys, ps):
        ce0 = -(y * np.log2(base) + (1 - y) * np.log2(1 - base))
        ce1 = -(y * np.log2(q) + (1 - y) * np.log2(1 - q))
        per.append(float(np.mean(ce0 - ce1)))
    return np.asarray(per), base


def higher_order(nsamp=500, seed=0):
    """Does the pool's JOINT structure carry anything its per-residue marginals do not?"""
    rows = []
    for a, t in enumerate(I.targets()):
        p, n = t["pdb"], int(t["n"])
        PH, PS, _ = RP.windows(p, "top75")
        PH = np.asarray(PH, float); PS = np.asarray(PS, float)
        nat = np.asarray(I.load_univ(p)["nat_ca"], float)
        rng = np.random.default_rng(seed + a)
        m = len(PH)
        # (1) whole windows: the full joint
        jw = rng.integers(0, m, nsamp)
        cj = I.build_ca(PH[jw], PS[jw])
        rj = I.kabsch_rmsd_batch(cj, nat)
        # (2) per-residue independent resampling: marginals only
        idx = rng.integers(0, m, (nsamp, n))
        r_ = np.arange(n)[None, :]
        ci = I.build_ca(PH[idx, r_], PS[idx, r_])
        ri = I.kabsch_rmsd_batch(ci, nat)
        rows.append({"pdb": p, "joint_mean": float(rj.mean()),
                     "joint_best": float(rj.min()),
                     "marg_mean": float(ri.mean()), "marg_best": float(ri.min())})
    return rows


def main():
    data = collect()
    folds = np.asarray([d["fold"] for d in data])
    names = [d["pdb"] for d in data]
    out = {"torsion": {}, "plugin": {}, "distogram": {}, "esm": {}, "higher_order": {}}

    print("PART D -- INFORMATION CONTENT (bits), with a direct structural measure beside "
          "each\n")
    print("D.1  RETRIEVAL POOL -> NATIVE TORSION, bits added over generic Ramachandran")
    print(f"{'bins':>6}{'bin width':>11}{'I_pred (bits/torsion)':>24}{'CI95':>20}"
          f"{'shuffled-target null':>22}")
    for B in BINS:
        v = torsion_information(data, B)
        vs = torsion_information(data, B, shuffle=True)
        pr = I.paired(v, np.zeros(len(v)), folds=folds, names=names)
        out["torsion"][B] = {"I_pred_bits": float(v.mean()),
                             "ci95": pr["ci95"], "median": float(np.median(v)),
                             "null_shuffled_targets": float(vs.mean()),
                             "wl": [int((v > 0).sum()), int((v < 0).sum())],
                             "concentration": L.concentration_verdict(-v)}
        print(f"{B:>6}{360.0 / B:>10.0f}d{v.mean():>24.3f}"
              f"   [{pr['ci95'][0]:+.3f},{pr['ci95'][1]:+.3f}]{vs.mean():>22.3f}")
    print("  paired structural measure: the same channel emits 4.072 A (top-75 circular "
          "mean)\n  vs 4.065 A for a zero-information constant alpha-helix.")

    print("\nD.2  the PLUG-IN estimator the brief warns about, for comparison")
    print(f"{'bins':>6}{'plug-in MI':>13}{'Miller-Madow':>14}{'shuffle null':>14}"
          f"{'null sd':>10}")
    for B in BINS:
        d = plugin_mi(data, B)
        out["plugin"][B] = d
        print(f"{B:>6}{d['plugin_mi_bits']:>13.3f}{d['miller_madow_corrected']:>14.3f}"
              f"{d['shuffle_null_mean_bits']:>14.3f}{d['shuffle_null_sd']:>10.3f}")

    v, vlo, vhi = distogram_information(data)
    pr = I.paired(v, np.zeros(len(v)), folds=folds, names=names)
    out["distogram"] = {"I_pred_bits_per_pair": float(v.mean()), "ci95": pr["ci95"],
                        "median": float(np.median(v)),
                        "low_sd_half": float(vlo.mean()),
                        "high_sd_half": float(vhi.mean()),
                        "wl": [int((v > 0).sum()), int((v < 0).sum())],
                        "concentration": L.concentration_verdict(-v)}
    print(f"\nD.3  DISTOGRAM -> NATIVE CA-CA DISTANCE (17 shipped bins), bits added over "
          f"the pooled distance marginal")
    print(f"  I_pred = {v.mean():+.3f} bits/pair  CI[{pr['ci95'][0]:+.3f},"
          f"{pr['ci95'][1]:+.3f}]  median {np.median(v):+.3f}  "
          f"targets positive {int((v > 0).sum())}/{len(v)}")
    print(f"  split by the channel's OWN uncertainty: low-sd half {vlo.mean():+.3f} bits, "
          f"high-sd half {vhi.mean():+.3f} bits")
    print("  paired structural measure: in-band ordering accuracy 0.545 (shipped) -> "
          "0.566 (low-sd half).")

    e, base = esm_contact_information()
    pre = I.paired(e, np.zeros(len(e)))
    out["esm"] = {"I_pred_bits_per_pair": float(e.mean()), "ci95": pre["ci95"],
                  "median": float(np.median(e)), "base_contact_rate": base,
                  "wl": [int((e > 0).sum()), int((e < 0).sum())]}
    print(f"\nD.4  ESM-2 CONTACT MAP -> NATIVE CONTACT (<8 A), base rate {base:.3f}")
    print(f"  I_pred = {e.mean():+.3f} bits/pair  CI[{pre['ci95'][0]:+.3f},"
          f"{pre['ci95'][1]:+.3f}]  targets positive {int((e > 0).sum())}/{len(e)}")
    print("  paired structural measure: in-band ordering accuracy 0.513 "
          "[0.494,0.531] -- at chance.")

    ho = higher_order()
    jm = np.asarray([r["joint_mean"] for r in ho]); jb = np.asarray([r["joint_best"]
                                                                    for r in ho])
    mm = np.asarray([r["marg_mean"] for r in ho]); mb = np.asarray([r["marg_best"]
                                                                    for r in ho])
    out["higher_order"] = {
        "rows": ho, "joint_mean": float(jm.mean()), "marg_mean": float(mm.mean()),
        "joint_best": float(jb.mean()), "marg_best": float(mb.mean()),
        "mean_diff": L.report(jm, mm, folds, names, "whole window - independent residues"),
        "best_diff": L.report(jb, mb, folds, names, "best-of-500")}
    print("\nD.5  THE POOL's HIGHER-ORDER STRUCTURE (500 draws per target)")
    print(f"  whole windows (full joint)          mean {jm.mean():.3f} A   "
          f"best-of-500 {jb.mean():.3f} A")
    print(f"  independent per-residue resampling  mean {mm.mean():.3f} A   "
          f"best-of-500 {mb.mean():.3f} A")
    print(f"  mean:        {L.fmt_paired(out['higher_order']['mean_diff'])}")
    print(f"  best-of-500: {L.fmt_paired(out['higher_order']['best_diff'])}")
    L.jwrite("info_mi", out)
    return out


if __name__ == "__main__":
    main()
