"""L2 / instrument diagnostic: how wide is the DEPOSITED NMR ensemble for each target?

ORACLE DIAGNOSTIC, NOT A METHOD.  This uses native coordinates and can never be an inference
-time signal.  The question it answers is about the BENCHMARK, not about a predictor:

  the manifest's reference is "deposited coordinates, MODEL 1".  For 117/126 targets the
  entry is an NMR ensemble of many conformers.  If the deposited models disagree with each
  other by d angstrom, then scoring against model 1 carries an irreducible reference term.

Reports, per target: n_models, mean pairwise CA-RMSD across deposited models, and RMSD of
model 1 to the ensemble medoid.
"""
import sys, os, json, urllib.request, collections
import numpy as np

sys.path.insert(0, r"C:\Users\abena\Protein-Folding-Algorithm")
import core.pipeline as P

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pdbcache")
os.makedirs(CACHE, exist_ok=True)


def fetch(pdb):
    p = os.path.join(CACHE, f"{pdb}.pdb")
    if os.path.exists(p) and os.path.getsize(p) > 200:
        return open(p, "r", errors="ignore").read()
    url = f"https://files.rcsb.org/download/{pdb}.pdb"
    with urllib.request.urlopen(url, timeout=90) as r:
        txt = r.read().decode("utf-8", "ignore")
    open(p, "w").write(txt)
    return txt


def models_ca(txt):
    """-> {chain: [ (model_idx, [(resseq, xyz), ...]) ]} for CA atoms, altloc ' '/'A'."""
    per = collections.defaultdict(lambda: collections.defaultdict(list))
    mi = 0
    for line in txt.splitlines():
        if line.startswith("MODEL"):
            try:
                mi = int(line[10:14])
            except ValueError:
                mi += 1
        elif line.startswith("ATOM") and line[12:16].strip() == "CA":
            alt = line[16]
            if alt not in (" ", "A"):
                continue
            ch = line[21]
            rs = line[22:27].strip()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            per[ch][mi].append((rs, xyz))
    return per


def kabsch_rmsd(A, B):
    A = A - A.mean(0); B = B - B.mean(0)
    U, S, Vt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return float(np.sqrt(((A @ R - B) ** 2).sum() / len(A)))


def main():
    targets = P.manifest("tuning126")
    meth = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                       "expmethod.json")))
    rows = []
    for t in targets:
        pdb = t.pdb.upper()
        n = len(t.ca)
        try:
            txt = fetch(pdb)
        except Exception as e:
            rows.append(dict(pdb=pdb, n=n, nmod=0, spread=np.nan, m1=np.nan,
                             method=meth.get(pdb, {}).get("method", "?"), err=str(e)[:40]))
            continue
        per = models_ca(txt)
        best = None
        for ch, mods in per.items():
            ks = sorted(mods)
            if not ks:
                continue
            L = len(mods[ks[0]])
            if L < n:
                continue
            # find the window of this chain matching model 1's deposited CA coords
            first = np.array([x for _, x in mods[ks[0]]], float)
            off, bd = None, 1e9
            for s in range(0, L - n + 1):
                d = np.abs(first[s:s + n] - t.ca).max()
                if d < bd:
                    bd, off = d, s
            if bd > 0.05:
                continue
            stack = []
            for k in ks:
                if len(mods[k]) != L:
                    continue
                stack.append(np.array([x for _, x in mods[k]], float)[off:off + n])
            if len(stack) >= 1 and (best is None or len(stack) > len(best)):
                best = stack
        if best is None or len(best) < 2:
            rows.append(dict(pdb=pdb, n=n, nmod=(0 if best is None else len(best)),
                             spread=np.nan, m1=np.nan,
                             method=meth.get(pdb, {}).get("method", "?"), err=""))
            continue
        X = np.array(best)
        M = len(X)
        D = np.zeros((M, M))
        for i in range(M):
            for j in range(i + 1, M):
                D[i, j] = D[j, i] = kabsch_rmsd(X[i], X[j])
        spread = float(D[np.triu_indices(M, 1)].mean())
        med = int(np.argmin(D.sum(1)))
        rows.append(dict(pdb=pdb, n=n, nmod=M, spread=spread, m1=float(D[0, med]),
                         method=meth.get(pdb, {}).get("method", "?"), err=""))

    ok = [r for r in rows if np.isfinite(r["spread"])]
    print(f"resolved ensembles for {len(ok)}/{len(rows)} targets\n")
    sp = np.array([r["spread"] for r in ok])
    m1 = np.array([r["m1"] for r in ok])
    print(f"mean pairwise CA-RMSD BETWEEN DEPOSITED MODELS : mean {sp.mean():.4f}  "
          f"median {np.median(sp):.4f}  sd {sp.std(ddof=1):.4f}  max {sp.max():.4f}")
    print(f"model 1 -> ensemble medoid                     : mean {m1.mean():.4f}  "
          f"median {np.median(m1):.4f}  max {m1.max():.4f}")
    print(f"targets with ensemble spread > 1.0 A : {(sp > 1.0).sum()}/{len(ok)}")
    print(f"targets with ensemble spread > 2.0 A : {(sp > 2.0).sum()}/{len(ok)}")
    print(f"targets with ensemble spread > 3.0 A : {(sp > 3.0).sum()}/{len(ok)}")
    print("\nwidest 15:")
    for r in sorted(ok, key=lambda r: -r["spread"])[:15]:
        print(f"  {r['pdb']}  n={r['n']:3d} models={r['nmod']:3d} "
              f"spread={r['spread']:7.3f}  m1->medoid={r['m1']:6.3f}  {r['method']}")
    json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "ens_spread.json"), "w"), indent=1)
    print("\nwrote ens_spread.json")


if __name__ == "__main__":
    main()
