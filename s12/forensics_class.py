"""forensics Part 2a: what do the near-native windows share, and is it deployable?

For every target: take the universe's 20 nearest windows (ORACLE selection) and ask what
their PARENTS have in common with the target that a deployable key could have seen:
  - structural class from the parent's own PDB header (lasso / fibril / membrane / ...),
    where the parent is a peptide (fragments have no peptide header)
  - sequence identity, BLOSUM rank, ESM-key rank, CA-SS agreement, hydrophobic pattern
  - the rank each of those keys would have given them
Control: the same statistics for 20 RANDOM universe windows and for the 20 windows the
shipped score ranks best.  Null: shuffled target-class labels.
"""
import os, sys, json, time, collections
os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.stats import mannwhitneyu, fisher_exact
from s12 import instrument as I
from s12 import forensics_lib as L

CLASSES = ["lasso", "fibril", "membrane", "cyclic", "bound", "designed", "cosolvent", "xray"]

def parent_headers():
    path = os.path.join(I.CACHE, "forensics_parent_headers.json")
    if os.path.exists(path):
        return json.load(path and open(path))
    P = L.parents(); out = {}
    for k in range(P["n_pep"]):
        out[str(P["pdb"][k])] = {c: bool(L.header_class(str(P["pdb"][k])).get(c, False)) for c in CLASSES}
    with open(path, "w") as fh:
        json.dump(out, fh)
    return out

def main():
    PH = parent_headers(); P = L.parents()
    tg = I.targets(); rows = []
    rng = np.random.default_rng(0)
    for k, t in enumerate(tg):
        pdb, n = t["pdb"], t["n"]
        u = I.load_univ(pdb); rr = u["rr"]
        pid, start = L.parent_map(pdb, u)
        cass = L.window_cass(pdb, u); esm = L.esm_key_cached(pdb, u)
        dg = I.distogram(pdb); i, j = I.pair_index(n)
        pool = I.pool_idx(u)
        D = I.pair_dists(u["W"][pool], i, j).astype(np.float32).astype(float)
        sc = I.shipped_score(dg, D)
        brank = np.empty(len(rr), int); brank[u["order"]] = np.arange(len(rr))
        erank = np.empty(len(rr), int); erank[np.argsort(-esm, kind="stable")] = np.arange(len(rr))
        tcls = L.header_class(pdb)
        near = np.argsort(rr, kind="stable")[:20]
        rand = rng.choice(len(rr), 20, replace=False)
        best_sc = pool[np.argsort(sc, kind="stable")[:20]]
        o_nss = L.ca_ss(u["nat_ca"])
        tcode = np.array([I.ALPHABET.index(c) for c in t["seq"]])
        def stats(idx, tag):
            pp = [str(P["pdb"][pid[w]]) for w in idx]
            ispep = np.array([bool(P["is_pep"][pid[w]]) for w in idx])
            shares = {}
            for c in CLASSES:
                if not tcls.get(c):
                    continue
                v = [PH.get(q, {}).get(c, False) for q, ip in zip(pp, ispep) if ip]
                shares[c] = float(np.mean(v)) if v else None
            return {f"{tag}_brank_med": float(np.median(brank[idx])), f"{tag}_erank_med": float(np.median(erank[idx])),
                    f"{tag}_in_pool": float((brank[idx] < I.K).mean()), f"{tag}_org": float(ispep.mean()),
                    f"{tag}_ident": float((u["S"][idx] == tcode[None, :]).mean()),
                    f"{tag}_ss_agree": float((cass[idx] == o_nss[None, :]).mean()),
                    f"{tag}_class_share": shares, f"{tag}_n_distinct_parents": len(set(pp)),
                    f"{tag}_rr": float(rr[idx].mean())}
        r = dict(pdb=pdb, fail18=pdb in I.FAIL18, n=n, tclass=[c for c in CLASSES if tcls.get(c)])
        r.update(stats(near, "near")); r.update(stats(rand, "rand")); r.update(stats(best_sc, "score"))
        r["near_parents"] = [(str(P["pdb"][pid[w]]), int(P["is_pep"][pid[w]]), float(rr[w]), int(brank[w]), int(erank[w])) for w in near[:8]]
        rows.append(r)
        print(f"[{k+1:3d}/126] {pdb}", flush=True)
    F = np.array([r["fail18"] for r in rows])
    agg = {}
    for m in ["near_brank_med", "near_erank_med", "near_in_pool", "near_org", "near_ident", "near_ss_agree", "near_n_distinct_parents",
              "rand_brank_med", "rand_erank_med", "rand_ident", "rand_ss_agree", "score_brank_med", "score_erank_med", "score_ident",
              "score_ss_agree", "score_rr", "near_rr"]:
        v = np.array([r[m] for r in rows], float)
        agg[m] = dict(all=float(v.mean()), fail18=float(v[F].mean()), other=float(v[~F].mean()), p=float(mannwhitneyu(v[F], v[~F]).pvalue))
    # class affinity: among targets carrying class c, does the near set over-represent peptide parents of class c?
    cls = {}
    base = {c: np.mean([PH[q][c] for q in PH]) for c in CLASSES}
    for c in CLASSES:
        got = [r for r in rows if c in r["tclass"] and r["near_class_share"].get(c) is not None]
        if len(got) < 3:
            continue
        near = np.array([r["near_class_share"][c] for r in got])
        rand = np.array([r["rand_class_share"][c] for r in got if r["rand_class_share"].get(c) is not None])
        cls[c] = dict(n_targets=len(got), near_share=float(near.mean()), rand_share=float(np.mean(rand)) if len(rand) else None,
                      library_base_rate=float(base[c]), enrichment=float(near.mean() / max(base[c], 1e-6)),
                      p=float(mannwhitneyu(near, rand).pvalue) if len(rand) else None)
    agg["class_affinity"] = cls; agg["library_base_rates"] = base
    I.write("forensics_class", dict(rows=rows, aggregate=agg))
    for m, v in agg.items():
        if m in ("class_affinity", "library_base_rates"):
            continue
        print(f"{m:28s} all {v['all']:9.3f}  F18 {v['fail18']:9.3f}  O108 {v['other']:9.3f}  p={v['p']:.3g}")
    print(json.dumps(cls, indent=1))
    print("base", {k: round(v, 3) for k, v in base.items()})

if __name__ == "__main__":
    main()
