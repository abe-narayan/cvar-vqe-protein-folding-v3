"""forensics: shared helpers.  Parent mapping of universe windows, per-window SS, keys,
the downstream chain (pool -> shipped score top-75 -> coordinate average -> project).

Nothing here touches benchmark60 or dev24.  `rr`/`nat_ca` are only read inside functions
whose names start with `oracle_` or in fields explicitly labelled as evaluation.
"""
from __future__ import annotations
import os, sys, json, re, time
import numpy as np
os.environ.setdefault("OMP_NUM_THREADS", "2")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I

CACHE = I.CACHE
SSMAP = {"H": 0, "E": 1, "C": 2}
SSINV = "HEC"
HYDRO = set("AVILMFWYC")

# ----------------------------------------------------------------------------- parents
_PARENTS = None
def parents():
    """Global parent table: peptides (787) then small fragments (6003)."""
    global _PARENTS
    if _PARENTS is None:
        import peptide_db as db, fragment_db as fdb
        peps = list(db.load()); frs = list(fdb.load(False))
        folds = db.folds(5)
        objs = peps + frs
        _PARENTS = dict(objs=objs, seq=np.array([q.seq for q in objs], dtype=object),
                        pdb=np.array([q.pdb for q in objs], dtype=object),
                        is_pep=np.array([True] * len(peps) + [False] * len(frs)),
                        fold=np.array([folds[q.seq] for q in peps] + [-1] * len(frs)),
                        n_pep=len(peps), folds=folds)
    return _PARENTS

def parent_map(pdb, u=None, verify=True):
    """(pid, start) for every window of the target's universe, reconstructed from the
    generation order and VERIFIED against the window's codes (exact) and CA."""
    path = os.path.join(CACHE, f"forensics_parents_{pdb}.npz")
    if os.path.exists(path):
        z = np.load(path); return z["pid"], z["start"]
    from core import data as cdata
    P = parents()
    if u is None:
        u = I.load_univ(pdb)
    fold, n, tseq = u["fold"], u["n"], u["seq"]
    pep_ids = [k for k in range(P["n_pep"]) if P["fold"][k] != fold and P["seq"][k] != tseq]
    fr = cdata.fold_fragments(fold, 5)
    frseq_to_id = {}
    for k in range(P["n_pep"], len(P["seq"])):
        frseq_to_id.setdefault((P["pdb"][k], P["seq"][k]), k)
    fr_ids = [frseq_to_id[(f.pdb, f.seq)] for f in fr]
    pid, start = [], []
    for k in pep_ids + fr_ids:
        m = len(P["seq"][k])
        if m < n:
            continue
        for s in range(m - n + 1):
            pid.append(k); start.append(s)
    pid = np.array(pid, np.int32); start = np.array(start, np.int16)
    assert len(pid) == len(u["W"]), (pdb, len(pid), len(u["W"]))
    if verify:
        codes = u["S"]
        expect = np.array([[I.ALPHABET.index(c) for c in P["seq"][pid[w]][start[w]:start[w] + n]]
                           for w in range(len(pid))], np.int8)
        assert (expect == codes).all(), (pdb, "full code mismatch")
        bad = 0
        for w in range(0, len(pid), max(1, len(pid) // 400)):
            ca = np.asarray(P["objs"][pid[w]].ca, float)[start[w]:start[w] + n]
            if np.abs(ca - u["W"][w]).max() > 1e-2:
                bad += 1
        assert bad == 0, (pdb, "parent CA mismatch", bad)
    np.savez_compressed(path, pid=pid, start=start)
    return pid, start

# ----------------------------------------------------------------------------- SS
def window_ss(pdb, u=None):
    """(nw, n) int8 H/E/C per window from the window's own torsions (I.ss_of).  Cached."""
    path = os.path.join(CACHE, f"forensics_ss_{pdb}.npz")
    if os.path.exists(path):
        return np.load(path)["ss"]
    if u is None:
        u = I.load_univ(pdb)
    out = np.zeros((len(u["W"]), u["n"]), np.int8)
    for w in range(len(u["W"])):
        out[w] = [SSMAP[c] for c in I.ss_of(u["PHI"][w], u["PSI"][w])]
    np.savez_compressed(path, ss=out)
    return out

def ca_ss(X):
    """CA-only SS (P-SEA-like): per residue H if d(i,i+3)~5.0-5.7 and d(i,i+4)~5.6-6.8 at
    i or i-1 (two consecutive turns), E if d(i,i+2)>6.3 and d(i,i+4)>11.0, else C.
    Works identically on nat_ca and on windows; used for every SS comparison."""
    X = np.asarray(X, float); n = len(X)
    d = lambda k: np.linalg.norm(X[k:] - X[:-k], axis=1) if k < n else np.zeros(0)
    d2, d3, d4 = d(2), d(3), d(4)
    turn = np.zeros(n, bool)
    for i in range(n - 4):
        if 4.9 <= d3[i] <= 5.7 and 5.6 <= d4[i] <= 6.8 and 4.9 <= d3[i + 1] <= 5.7:
            turn[i] = True
    ss = np.full(n, 2, np.int8)
    for i in range(n - 4):
        if turn[i]:
            ss[i:i + 5] = 0
    for i in range(n - 4):
        if d2[i] > 6.3 and d2[i + 2] > 6.3 and d4[i] > 11.0 and (ss[i:i + 5] == 2).all():
            ss[i:i + 5] = 1
    return ss

def ca_ss_batch(W):
    return np.stack([ca_ss(w) for w in W])

def window_cass(pdb, u=None):
    """(nw, n) int8 CA-geometry SS per window.  Cached."""
    path = os.path.join(CACHE, f"forensics_cass_{pdb}.npz")
    if os.path.exists(path):
        return np.load(path)["ss"]
    if u is None:
        u = I.load_univ(pdb)
    out = ca_ss_batch(u["W"]).astype(np.int8)
    np.savez_compressed(path, ss=out)
    return out

def esm_key_cached(pdb, u=None):
    path = os.path.join(CACHE, f"forensics_esm_{pdb}.npz")
    if os.path.exists(path):
        return np.load(path)["key"]
    if u is None:
        u = I.load_univ(pdb)
    pid, start = parent_map(pdb, u)
    key = esm_key(pdb, u, pid, start)
    np.savez_compressed(path, key=key)
    return key

def ss_str(a):
    return "".join(SSINV[int(x)] for x in np.asarray(a).ravel())

def oracle_native_ss(pdb):
    """Native SS from the target's own deposited torsions (LABEL ONLY)."""
    import peptide_db as db
    p = db.by_pdb(pdb)
    return np.array([SSMAP[c] for c in I.ss_of(p.phi, p.psi)], np.int8), p

def distogram_ss(dg, n, h4=6.5, e2=6.4):
    """Deployable SS-like class per residue from the shipped distogram's expected distances.
    helix if E[d(i,i+4)] < h4 and E[d(i+1,i+4)] < 6.2 (assigned to i..i+4); strand if
    E[d(i,i+2)] > e2 on non-helical residues."""
    i, j = dg["i"], dg["j"]; ex = dg["expected"]
    d = {(int(a), int(b)): float(e) for a, b, e in zip(i, j, ex)}
    ss = np.full(n, 2, np.int8)
    for a in range(n - 4):
        if d.get((a, a + 4), 99) < h4 and d.get((a + 1, a + 4), 99) < 6.2:
            ss[a:a + 5] = 0
    for a in range(n - 2):
        if ss[a] == 2 and ss[a + 2] == 2 and d.get((a, a + 2), 0) > e2:
            ss[a:a + 3] = np.where(ss[a:a + 3] == 2, 1, ss[a:a + 3])
    return ss

# ----------------------------------------------------------------------------- keys
def esm_key(pdb, u, pid, start):
    """Mean per-position cosine (pca128) between target residues and window residues."""
    from s12 import esm_bank
    bank = esm_bank.load(); P = parents()
    n = u["n"]
    T = bank[u["seq"]][1]; T = T / np.linalg.norm(T, axis=1, keepdims=True)
    key = np.zeros(len(pid), np.float32)
    cache = {}
    for w in range(len(pid)):
        k = int(pid[w])
        if k not in cache:
            E = bank[P["seq"][k]][1]; cache[k] = E / np.linalg.norm(E, axis=1, keepdims=True)
        E = cache[k][start[w]:start[w] + n]
        key[w] = (E * T).sum(1).mean()
    return key

def hydro_pattern(seq):
    return np.array([c in HYDRO for c in seq], bool)

def rank_key(key, k=I.K):
    """Stable argsort of -key (ties resolved by universe position, like the shipped pool)."""
    return np.argsort(-np.asarray(key, float), kind="stable")[:k]

# ----------------------------------------------------------------------------- chain
def chain(u, pool, dg, seq, fold, do_project=True, m=I.M):
    """pool (indices into the universe) -> shipped-score top-m -> coordinate average -> project.
    rr/nat_ca are used ONLY for the evaluation fields."""
    n = u["n"]; i, j = I.pair_index(n)
    pool = np.asarray(pool, int)
    Wp = u["W"][pool]
    D = I.pair_dists(Wp, i, j).astype(np.float32).astype(float)
    sc = I.shipped_score(dg, D)
    top = np.argsort(sc, kind="stable")[:m]
    C, b = I.coordinate_average(Wp[top])
    out = {"pool_size": int(len(pool)), "top_idx_univ": pool[top].tolist()}
    rr = u["rr"]
    out["pool_best"] = float(rr[pool].min()); out["pool_mean"] = float(rr[pool].mean())
    out["argmin"] = float(rr[pool[int(np.argmin(sc))]])
    out["top_best"] = float(rr[pool[top]].min()); out["top_mean"] = float(rr[pool[top]].mean())
    out["rmsd_avg"] = I.ca_rmsd(C, u["nat_ca"])
    if do_project:
        pj = I.project(C, seq, fold)
        out["rmsd_fit"] = I.ca_rmsd(pj["fit_ca"], u["nat_ca"]); out["rmsd_arm"] = I.ca_rmsd(pj["ca"], u["nat_ca"])
    return out

# ----------------------------------------------------------------------------- headers
_HDR_RULES = [
    ("lasso", r"LASSO"),
    ("cyclic", r"CYCLIC|CYCLO|HEAD-TO-TAIL|MACROCYCL"),
    ("fibril", r"FIBRIL|AMYLOID|STERIC ZIPPER"),
    ("xray", r"^EXPDTA\s+(X-RAY|ELECTRON CRYST)"),
    ("membrane", r"MICELLE|\bSDS\b|\bDPC\b|BICELLE|MEMBRANE|LIPID|DHPC|DLPC|TRANSMEMBRANE"),
    ("cosolvent", r"\bTFE\b|TRIFLUOROETHANOL|DMSO|DIMETHYL\s*SULFOXIDE|HFIP|METHANOL"),
    ("bound", r"\bBOUND\b|IN COMPLEX|COMPLEX WITH|BINDING TO|TRANSFERRED NOE|EXCHANGE TRANSFER"),
    ("designed", r"DE NOVO|DESIGNED|SYNTHETIC"),
]
def header_class(pdb):
    """Metadata-only classification from the PDB header (HEADER/TITLE/KEYWDS/EXPDTA/REMARK 210/SSBOND/LINK)."""
    path = None
    for d in ("pdbs_ext", "pdbs"):
        q = os.path.join(ROOT, d, f"{pdb}.pdb")
        if os.path.exists(q):
            path = q; break
    if path is None:
        return {"found": False}
    hdr, ssbond, link_iso, link_cyc, nmodel = [], 0, 0, 0, 0
    with open(path, errors="ignore") as fh:
        for line in fh:
            if line.startswith(("HEADER", "TITLE", "KEYWDS", "EXPDTA", "COMPND", "REMARK 210")):
                hdr.append(line.rstrip())
            elif line.startswith("SSBOND"):
                ssbond += 1
            elif line.startswith("LINK"):
                a1, a2 = line[12:16].strip(), line[42:46].strip()
                r1, r2 = line[17:20].strip(), line[47:50].strip()
                caps = {"NH2", "ACE", "NME", "HOH"}
                if r1 in caps or r2 in caps:
                    continue
                if (a1 == "N" and a2 in ("CG", "CD") and r2 in ("ASP", "GLU")) or \
                   (a2 == "N" and a1 in ("CG", "CD") and r1 in ("ASP", "GLU")):
                    link_iso += 1
                elif {a1, a2} == {"N", "C"}:
                    link_cyc += 1
            elif line.startswith("MODEL"):
                nmodel += 1
    text = "\n".join(hdr)
    out = {"found": True, "ssbond": ssbond, "link_isopeptide": link_iso, "link_cyclic": link_cyc, "n_models": nmodel,
           "title": " ".join(l[10:].strip() for l in hdr if l.startswith("TITLE"))[:160],
           "expdta": " ".join(l[10:].strip() for l in hdr if l.startswith("EXPDTA"))[:60]}
    for name, rx in _HDR_RULES:
        out[name] = bool(re.search(rx, text, flags=re.M))
    out["lasso"] = out["lasso"] or link_iso > 0
    out["cyclic"] = out["cyclic"] or link_cyc > 0
    out["covalent"] = out["lasso"] or out["cyclic"] or ssbond > 0
    out["nonaqueous"] = out["membrane"] or out["cosolvent"] or out["xray"]
    out["any_flag"] = out["covalent"] or out["nonaqueous"] or out["fibril"] or out["bound"]
    return out

def rg(X):
    X = np.asarray(X, float); return float(np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean()))

def cluster_count(P, cut=2.5):
    """Number of clusters (average linkage at `cut` A) of a pairwise RMSD matrix."""
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    P = np.asarray(P, float); P = 0.5 * (P + P.T); np.fill_diagonal(P, 0.0)
    if len(P) < 2:
        return 1, np.zeros(len(P), int)
    Z = linkage(squareform(P, checks=False), "average")
    lab = fcluster(Z, cut, criterion="distance")
    return int(lab.max()), lab
