"""EXPERIMENT 3 -- MULTI-HYPOTHESIS OUTPUT.

Cluster the shipped top-75 into k structural modes, emit a per-mode coordinate average, and
ask (a) how much ORACLE headroom best-of-k has, (b) whether ANY native-free rule for picking
or weighting the modes reaches any of it, and (c) whether a MODE-AWARE average (equal weight
per mode, i.e. de-biasing the average for mode population) beats the global average even
when the mode cannot be picked.

Clusterings: agglomerative average / complete / ward on the candidate-vs-candidate CA-RMSD
matrix, and k-medoids, for k = 2..5.

Native-free mode-picking rules (13): largest, tightest, best mean score, best score of the
mode average, best risk of the mode's mean pair distances, contains the argmin, closest to
the global average, best BLOSUM sim, most peptide-derived (org), best legacy energy, best
mean rank, smallest |rg - global rg|, and a soft mixture over modes.
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import agg_features as AF
from s12 import agg_meta
from s12 import agg_subset as SS

CLOUD = os.path.join(ROOT, "s12", "cache", "agg_modes")
os.makedirs(CLOUD, exist_ok=True)


def cluster(P, k, method):
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    if method == "kmed":
        rng = np.random.default_rng(0)
        m = len(P)
        cen = [int(np.argmin(P.mean(1)))]
        while len(cen) < k:
            d = P[cen].min(0)
            cen.append(int(np.argmax(d)))
        for _ in range(30):
            lab = np.argmin(P[cen], 0)
            new = []
            for c in range(k):
                ix = np.where(lab == c)[0]
                if len(ix) == 0:
                    new.append(cen[c]); continue
                new.append(int(ix[np.argmin(P[np.ix_(ix, ix)].mean(1))]))
            if new == cen:
                break
            cen = new
        return np.argmin(P[cen], 0)
    Z = linkage(squareform(P, checks=False), method=method)
    return fcluster(Z, t=k, criterion="maxclust") - 1


def mode_avgs(A, P, lab):
    out = []
    for c in range(lab.max() + 1):
        ix = np.where(lab == c)[0]
        if len(ix) == 0:
            continue
        out.append((ix, A[ix].mean(0)))
    return out


def run(out="agg_modes"):
    tg = I.targets()
    rows = {}
    t0 = time.time()
    for q, t in enumerate(tg):
        ins = SS.Inst(t["pdb"])
        mt = agg_meta.load(t["pdb"])
        o = AF.load(t["pdb"])
        idx75 = o["idx"].astype(int)
        A, P, nat, sc = ins.A, ins.P, ins.nat, ins.sc
        Cg = A.mean(0)
        rg_g = np.sqrt(((Cg - Cg.mean(0)) ** 2).sum(-1).mean())
        r = {"avg75": I.ca_rmsd(Cg, nat), "n": t["n"]}
        for method in ("average", "ward", "kmed"):
            for k in (2, 3, 4, 5):
                lab = cluster(P, k, method)
                ma = mode_avgs(A, P, lab)
                Cs = np.stack([c for _, c in ma])
                rms = I.kabsch_rmsd_batch(Cs, nat)
                tag = f"{method}{k}"
                r[f"{tag}_ORC_best"] = float(rms.min())
                r[f"{tag}_ORC_worst"] = float(rms.max())
                r[f"{tag}_nmodes"] = len(ma)
                r[f"{tag}_sizes"] = [int(len(ix)) for ix, _ in ma]
                # ---- mode-aware average (equal weight per mode) -- needs no mode pick
                r[f"{tag}_modeavg"] = I.ca_rmsd(Cs.mean(0), nat)
                # ---- native-free picking rules
                feats = {}
                feats["largest"] = np.array([-len(ix) for ix, _ in ma], float)
                feats["tightest"] = np.array([P[np.ix_(ix, ix)].mean() for ix, _ in ma])
                feats["meanscore"] = np.array([sc[ix].mean() for ix, _ in ma])
                feats["minscore"] = np.array([sc[ix].min() for ix, _ in ma])
                feats["riskavg"] = np.array([ins.risk_of(c) for _, c in ma])
                feats["riskdist"] = np.array([ins.risk_dist(ix) for ix, _ in ma])
                feats["hasargmin"] = np.array([-float(int(np.argmin(sc)) in set(ix.tolist())) for ix, _ in ma])
                feats["nearglobal"] = np.array([I.ca_rmsd(c, Cg) for _, c in ma])
                feats["blosum"] = np.array([-mt["sim"][idx75][ix].mean() for ix, _ in ma])
                feats["org"] = np.array([-mt["org"][idx75][ix].mean() for ix, _ in ma])
                feats["rgclose"] = np.array([abs(np.sqrt(((c - c.mean(0)) ** 2).sum(-1).mean()) - rg_g) for _, c in ma])
                feats["meanrank"] = np.array([np.argsort(np.argsort(sc))[ix].mean() for ix, _ in ma])
                feats["spread"] = np.array([np.mean([I.ca_rmsd(A[i], c) for i in ix]) for ix, c in ma])
                for nmn, fv in feats.items():
                    # average over the tied argmin set (avoids the np.argmin tie-break leak)
                    lo = np.where(fv <= fv.min() + 1e-12)[0]
                    r[f"{tag}_pick_{nmn}"] = float(rms[lo].mean())
                # soft mixture over modes by riskdist
                z = (feats["riskdist"] - feats["riskdist"].mean()) / (feats["riskdist"].std() + 1e-9)
                for T in (0.5, 1.0):
                    w = np.exp(-z / T); w /= w.sum()
                    r[f"{tag}_soft{T}"] = I.ca_rmsd((Cs * w[:, None, None]).sum(0), nat)
        rows[t["pdb"]] = r
        if q % 10 == 0:
            print(q, t["pdb"], f"{time.time()-t0:.0f}s", flush=True)
    I.write(out, rows)
    print("done", time.time() - t0)


if __name__ == "__main__":
    run()
