"""LANE L -- `long40`, a SEPARATE 40-60 residue instrument.  Built to PREREG S32 L, section 3.

`long40` is never merged with, re-split against, or substituted for `tuning126`, and no
S32 headline about the 3.2105 A endpoint may be computed on it (contract rule 11).  It has
its OWN fold structure, frozen before the first arm runs.

ADMISSION, in the fixed order the pre-registration states:
  1. `prots/` chain, first model, first chain, 40 <= n <= 60
  2. contiguous CA (3.5-4.1 A), finite interior torsions, no X          [L1 census gate]
  3. single-chain deposition -- a monomer, so the deposited conformation is not held in
     place by an interface the prediction cannot see
  4. LEAKAGE: the target's own PDB contributes nothing to `peptide_db` or `fragment_db`,
     and no bank member appears VERBATIM inside the target covering >= 0.25 of it.

     THE FILTER WAS CHANGED AFTER MEASURING ITS NULL, and the first choice was wrong.
     The pre-registration named `norm="shorter"` identity >= 0.4.  Run, it rejected 140 of
     170 monomers and left ZERO targets.  `s32_L_leaknull` then measured what the project's
     memory said to measure (`containment-threshold-is-at-the-null`): on 60 targets x 3
     composition-preserving shuffles each, the `shorter` statistic rejects **100% of REAL
     and 100% of SHUFFLED sequences at every threshold up to 0.9**.  It is exactly at the
     null and carries no information about leakage: normalised by a 9-residue fragment,
     4 chance matches already score 0.44.

     The verbatim-substring statistic separates cleanly instead -- real mean 0.330,
     shuffled mean **0.000**, shuffled max 0.000 -- so it is the statistic used.  This is a
     CORRECTION MADE ON A NULL MEASUREMENT, not on an outcome: no RMSD was computed on any
     target at the time it was made, and the threshold is set at the point the null allows
     (any positive threshold is infinitely above a shuffled rate of zero), not at the point
     that maximises survivors.
  5. REDUNDANCY: single-linkage cluster survivors at identity >= 0.4 (longer-normalised,
     the pinned convention, both sequences comparable in length here); keep one
     representative per cluster, the LOWEST PDB CODE -- a rule fixed in the
     pre-registration before any RMSD was known, so it cannot be chosen on outcome.

FOLDS: 5, assigned by greedy balanced packing of whole clusters, largest first, so no
cluster is split across folds.  Deterministic, seeded, written once.

    python -m s32.s32_L_instrument build
    python -m s32.s32_L_instrument bank --pdb 1ABC
"""
from __future__ import annotations

import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "s32", "results")
os.makedirs(RESULTS, exist_ok=True)

LO, HI = 40, 60
#: Verbatim-substring coverage of the target by any bank member.  Shuffled controls score
#: 0.000 at this statistic (n = 180 shuffles), so any positive threshold is above the null;
#: 0.25 keeps it strict without being arbitrary -- it excludes a target the moment a
#: quarter of it is a sequence the distance prior has already been trained on.
LEAK_VERBATIM = 0.25
#: Protein-to-protein redundancy: both sequences are 40-60 residues, so the two
#: normalisations coincide and the pinned `longer` convention is used unmodified.
CLUSTER_ID = 0.4
#: Retained only to exclude a target's own chain from its candidate bank.
LEAK_ID = 0.4
N_FOLDS = 5
MANIFEST = os.path.join(RESULTS, "long40_manifest.json")
FOLDS = os.path.join(RESULTS, "long40_folds.json")


def _nchains(path):
    chains = set()
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM") and line[12:16].strip() == "CA" \
                    and line[16] in (" ", "A"):
                chains.add(line[21])
    return len(chains)


def build(verbose=True):
    from core import data as cdata
    import peptide_db as pdb
    import fragment_db as fdb

    cen = json.load(open(os.path.join(RESULTS, "L1_census.json")))
    cand = [r for r in cen["rows"] if LO <= r["n"] <= HI]
    cand.sort(key=lambda r: r["pdb"])
    step = {"in_band": len(cand)}

    # --- gate 3: monomer
    keep = []
    for r in cand:
        p = os.path.join(BASE, "prots", r["pdb"] + ".pdb")
        if not os.path.exists(p):
            g = glob.glob(os.path.join(BASE, "prots", r["pdb"].lower() + ".pdb"))
            if not g:
                continue
            p = g[0]
        if _nchains(p) == 1:
            r = dict(r); r["path"] = p
            keep.append(r)
    cand = keep
    step["monomer"] = len(cand)

    # --- gate 4: leakage against the two banks the distogram fold models train on
    peps = list(pdb.load())
    frags = list(fdb.load())
    frag_src = {f.pdb.split("_")[0].upper() for f in frags}
    pep_src = {p.pdb.split("_")[0].upper() for p in peps}
    bank_set = {p.seq for p in peps} | {f.seq for f in frags}

    keep, why = [], {"own-pdb-in-bank": 0, "verbatim-bank-member": 0}
    for r in cand:
        if r["pdb"].upper() in frag_src or r["pdb"].upper() in pep_src:
            why["own-pdb-in-bank"] += 1
            continue
        vb = max((len(b) for b in bank_set if b in r["seq"]), default=0) / len(r["seq"])
        if vb >= LEAK_VERBATIM:
            why["verbatim-bank-member"] += 1
            continue
        r = dict(r); r["verbatim_bank_coverage"] = float(vb)
        keep.append(r)
    cand = keep
    step["leakage_survivors"] = len(cand)
    step["leakage_rejected"] = why
    if verbose:
        print(f"after leakage: {len(cand)}  ({why})", flush=True)

    # --- gate 5: redundancy, single linkage at CLUSTER_ID, lowest PDB code kept
    seqs = [r["seq"] for r in cand]
    n = len(cand)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a

    for a in range(n):
        idv = cdata.identity_many(seqs[a], seqs[a + 1:])
        for off, v in enumerate(idv):
            if v >= CLUSTER_ID:
                ra, rb = find(a), find(a + 1 + off)
                if ra != rb:
                    parent[ra] = rb
    groups = {}
    for a in range(n):
        groups.setdefault(find(a), []).append(a)
    reps = []
    for g in groups.values():
        g.sort(key=lambda a: cand[a]["pdb"])
        r = dict(cand[g[0]]); r["cluster_size"] = len(g)
        r["cluster_members"] = [cand[a]["pdb"] for a in g]
        reps.append(r)
    reps.sort(key=lambda r: r["pdb"])
    step["n_clusters"] = len(reps)

    # --- folds: greedy balanced packing of whole clusters, largest first
    order = sorted(range(len(reps)), key=lambda a: (-reps[a]["cluster_size"],
                                                    reps[a]["pdb"]))
    load = [0] * N_FOLDS
    for a in order:
        f = int(np.argmin(load))
        reps[a]["fold"] = f
        load[f] += 1
    out = dict(name="long40", band=[LO, HI], leak_identity=LEAK_ID,
               cluster_identity=CLUSTER_ID, n_folds=N_FOLDS,
               fold_sizes=load, funnel=step, n=len(reps),
               lengths=dict(mean=float(np.mean([r["n"] for r in reps])),
                            min=int(min(r["n"] for r in reps)),
                            max=int(max(r["n"] for r in reps))),
               targets=reps)
    with open(MANIFEST, "w") as fh:
        json.dump(out, fh, indent=1)
    with open(FOLDS, "w") as fh:
        json.dump({r["pdb"]: int(r["fold"]) for r in reps}, fh, indent=1)
    if verbose:
        print(json.dumps({k: v for k, v in out.items() if k != "targets"}, indent=1))
        print("wrote", MANIFEST, "and", FOLDS, flush=True)
    return out


def targets():
    return json.load(open(MANIFEST))["targets"]


# ------------------------------------------------------------------ the candidate bank
def bank(pdb, K=500, verbose=False, exclude_id=LEAK_ID):
    """Top-K BLOSUM62 length-n windows of `prots/`, leakage-filtered, STABLE argsort.

    Excludes the target's own chain and every chain reaching `exclude_id` identity to it.
    The window universe at n = 45 is 1.7 M, far too large to cache per target the way
    `s8/generate_univ` does at n = 13, so the scan is STREAMING: BLOSUM sums are
    accumulated per protein and only the running top-K coordinates are held.
    """
    from core import data as cdata
    from core import geometry as geo
    import s7.audit as audit

    tg = {t["pdb"]: t for t in targets()}[pdb]
    n = int(tg["n"]); tseq = tg["seq"]
    code = audit.encode(tseq)
    paths = sorted(glob.glob(os.path.join(BASE, "prots", "*.pdb")))
    best_sim, best_W, best_PH, best_PS, best_src = [], [], [], [], []
    t0 = time.time()
    for k, p in enumerate(paths):
        src = os.path.basename(p)[:-4].upper()
        if src == pdb.upper():
            continue
        try:
            seq, coords, phi, psi = geo.native_coords_from_pdb(p)
        except Exception:
            continue
        m = len(seq)
        if m < n or len(coords["CA"]) != m:
            continue
        if float(cdata.identity_many(tseq, [seq], norm="shorter")[0]) >= exclude_id:
            continue
        ca = np.asarray(coords["CA"], float)
        step = np.linalg.norm(np.diff(ca, axis=0), axis=1)
        ok = (step > 3.5) & (step < 4.1)
        phi = np.asarray(phi, float); psi = np.asarray(psi, float)
        starts = [s for s in range(1, m - n)
                  if ok[s:s + n - 1].all() and "X" not in seq[s:s + n]
                  and np.all(np.isfinite(phi[s + 1:s + n]))
                  and np.all(np.isfinite(psi[s:s + n - 1]))]
        if not starts:
            continue
        S = np.stack([audit.encode(seq[s:s + n]) for s in starts])
        sim = audit.B62[S, code[None, :]].sum(1)
        for a, s in enumerate(starts):
            best_sim.append(float(sim[a])); best_src.append(f"{src}_{s}")
            best_W.append(ca[s:s + n]); best_PH.append(phi[s:s + n])
            best_PS.append(psi[s:s + n])
        if len(best_sim) > 40 * K:                      # periodic prune, order-preserving
            idx = np.argsort(-np.asarray(best_sim), kind="stable")[:4 * K]
            idx = np.sort(idx)
            best_sim = [best_sim[a] for a in idx]; best_src = [best_src[a] for a in idx]
            best_W = [best_W[a] for a in idx]; best_PH = [best_PH[a] for a in idx]
            best_PS = [best_PS[a] for a in idx]
        if verbose and (k + 1) % 4000 == 0:
            print(f"    {k+1}/{len(paths)} scanned, {len(best_sim)} held "
                  f"({time.time()-t0:.0f}s)", flush=True)
    idx = np.argsort(-np.asarray(best_sim), kind="stable")[:K]
    return dict(pdb=pdb, n=n, seq=tseq,
                W=np.stack([best_W[a] for a in idx]).astype(np.float32),
                PHI=np.stack([best_PH[a] for a in idx]).astype(np.float32),
                PSI=np.stack([best_PS[a] for a in idx]).astype(np.float32),
                sim=np.asarray([best_sim[a] for a in idx], np.float32),
                src=[best_src[a] for a in idx], secs=time.time() - t0)


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "build"
    if stage == "build":
        build()
    elif stage == "bank":
        b = bank(sys.argv[sys.argv.index("--pdb") + 1], verbose=True)
        print(f"{b['pdb']} n={b['n']} K={len(b['sim'])} in {b['secs']:.0f}s")
