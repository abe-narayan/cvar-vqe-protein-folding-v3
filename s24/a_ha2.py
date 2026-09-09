"""SPRINT 24 / WORKSTREAM A / H-A2 -- is the SHIPPED fragment filter leaky, and how much?

Settled from source and from training-set identifiers only.  No benchmark result file and
no benchmark structure is opened; benchmark SEQUENCES and IDs are read solely to test
whether a training chain is related to one.

THE SOURCE, read rather than assumed
------------------------------------
`distogram._fold_fragments(fold, n_folds, threshold=0.6)`:

    held = [p.seq for p in pdb.load() if pdb.folds(n_folds)[p.seq] == fold]
    hk   = [pdb._kmers(s) for s in held]                        # 3-mer sets
    for f in fdb.load():
        fk = pdb._kmers(f.seq)
        for s, k in zip(held, hk):
            if fk and k and not (fk & k):     # <-- PREFILTER (A)
                continue
            if pdb.identity(s, f.seq) >= threshold:   # <-- normalised by the LONGER (B)
                leak = True

and `core/pipeline.fold_model(fold)` -> `train_fold(fold, True, 5, fragments=True)` trains
on out-of-fold peptides PLUS exactly this fragment set.  `core/bench.py:207` and
`core/pipeline.py:1315` both set `fold = db.folds(cfg.n_folds)[target.seq]`.

So the target's OWN fold is the held-out fold at inference, and the fragment filter IS run
against the target's own sequence.  Two residual holes remain, and they are what this
module sizes:

  (A) the 3-mer prefilter.  `peptide_db.clusters` removed exactly this prefilter as
      UNSOUND ("a pair can reach 0.6 identity with no shared 3-mer at all") after it let a
      benchmark target train on a 0.70-identity homolog.  `_fold_fragments` still has it.
  (B) identity normalised by the LONGER sequence.  A 9-mer target sitting VERBATIM inside a
      20-mer fragment scores 9/20 = 0.45 and passes a 0.60 threshold untouched.  This is
      the same defect memory records for `peptide_db` (4/126 tuning targets carry a verbatim
      copy in their own distogram's training set); fragments run to 20 aa against 9-16 aa
      targets, so the exposure is structurally larger here.

Usage: python -m s24.a_ha2
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import peptide_db as pdb                                                # noqa: E402
import fragment_db as fdb                                              # noqa: E402
import distogram as dgm                                                # noqa: E402
import core.data as cdata                                              # noqa: E402
from s7 import debias                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
T = pdb.IDENTITY_THRESHOLD


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True, default=str)
    os.replace(tmp, path)


def dep(pid):
    return str(pid).split("_")[0].upper()[:4]


def main():
    t0 = time.time()
    folds = pdb.folds(5)
    peps = list(pdb.load())
    frags = list(fdb.load())
    tun = debias.tuning_targets()
    dev = cdata.dev_set(24)
    bench = cdata.benchmark()
    sets = {"tuning126": tun, "dev24": dev, "bench60": bench}

    # the shipped fragment training set for each fold
    ffold = {f: list(dgm._fold_fragments(f, 5, T)) for f in range(5)}
    print("shipped fragment training set per fold: "
          + ", ".join(f"{f}:{len(ffold[f])}" for f in range(5))
          + f"  (of {len(frags)})", flush=True)

    # the shipped PEPTIDE training set for each fold
    pfold = {f: [p for p in peps if folds[p.seq] != f] for f in range(5)}

    out = {"n_fragments": len(frags), "n_peptides": len(peps),
           "frag_per_fold": {f: len(ffold[f]) for f in range(5)},
           "pep_per_fold": {f: len(pfold[f]) for f in range(5)},
           "sets": {}}

    for name, tgts in sets.items():
        rec = []
        for p in tgts:
            f = folds[p.seq]
            tr_f, tr_p = ffold[f], pfold[f]
            hit = {"pdb": p.pdb, "n": p.n, "fold": int(f),
                   "frag_verbatim": [], "frag_homolog_missed_by_kmer": [],
                   "frag_deposit": [], "pep_verbatim": [], "pep_deposit": []}
            tk = pdb._kmers(p.seq)
            ct = pdb._composition(p.seq)
            for q in tr_f:
                if p.seq in q.seq or q.seq in p.seq:     # (B) verbatim, EITHER direction
                    hit["frag_verbatim"].append(
                        [q.pdb, q.seq, round(pdb.identity(p.seq, q.seq), 3)])
                else:
                    qk = pdb._kmers(q.seq)
                    if tk and qk and not (tk & qk):      # (A) prefilter would have skipped
                        if pdb.max_possible_identity(p.seq, q.seq, ct) >= T \
                                and pdb.identity(p.seq, q.seq) >= T:
                            hit["frag_homolog_missed_by_kmer"].append(
                                [q.pdb, q.seq, round(pdb.identity(p.seq, q.seq), 3)])
                if dep(q.pdb) == dep(p.pdb):
                    hit["frag_deposit"].append([q.pdb, q.seq])
            for q in tr_p:
                if p.seq in q.seq or q.seq in p.seq:
                    hit["pep_verbatim"].append(
                        [q.pdb, q.seq, round(pdb.identity(p.seq, q.seq), 3)])
                if dep(q.pdb) == dep(p.pdb):
                    hit["pep_deposit"].append([q.pdb, q.seq])
            hit["any"] = bool(hit["frag_verbatim"] or hit["frag_homolog_missed_by_kmer"]
                              or hit["frag_deposit"] or hit["pep_verbatim"]
                              or hit["pep_deposit"])
            rec.append(hit)
        summ = {k: int(sum(bool(r[k]) for r in rec)) for k in
                ("frag_verbatim", "frag_homolog_missed_by_kmer", "frag_deposit",
                 "pep_verbatim", "pep_deposit", "any")}
        summ["n_targets"] = len(rec)
        out["sets"][name] = {"summary": summ,
                             "affected": [r for r in rec if r["any"]],
                             "complete": len(rec) == len(tgts)}
        print(f"{name:10s} n={len(rec):3d}  " + json.dumps(summ), flush=True)

    out["complete"] = all(v["complete"] for v in out["sets"].values())
    write_json(os.path.join(RESULTS, "a_ha2_fragment_leak.json"), out)
    print(f"{time.time()-t0:.0f}s -> s24/results/a_ha2_fragment_leak.json")


if __name__ == "__main__":
    main()
