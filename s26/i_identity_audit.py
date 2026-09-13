#!/usr/bin/env python
"""s26/i_identity_audit.py -- what the corrected identity normalisation would do to the folds,
computed IN MEMORY.  Sprint 26, lane I, declared defect 6(a).

The defect (`core/data.py:182`, `peptide_db.py:119`): Needleman-Wunsch identity is normalised
by the LONGER sequence, so a 13-mer sitting verbatim inside a 25-mer scores 13/25 = 0.520,
passes the 0.6 clustering threshold, and can land in a different fold from its own copy.  The
corrected measure is the same alignment normalised by the SHORTER sequence
(`core.data.containment`, exposed as `core.data.identity(..., norm="shorter")` since S26) with
an explicit verbatim-substring test in front of it.

What this script never does: it never calls `core.data.clusters()` / `folds()` or their
`peptide_db` twins, which are write-on-first-use caches, and it never opens
`peptide_folds.json` or `peptide_clusters.json` for writing.  The clustering logic of
`core.data.clusters` is copied here and run twice in memory -- once under the pinned convention
(and checked against the pinned files, so the copy is proven faithful) and once under the
corrected one.  Both pinned files are sha256-hashed before and after and compared to lane E's
`s26/results/pinned_hashes.json`.  The sealed benchmark is touched only through
`core.backend("data").benchmark()` for its SEQUENCES; no benchmark name and no RMSD is read,
printed or saved.

Output: `s26/results/i_identity_audit.json` and a printed summary.

    python s26/jobrun.py --agent I --tag CPU --name i_identity_audit --est-ram 0.5 \
        -- python s26/i_identity_audit.py
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import core                                                              # noqa: E402
from core import data as D                                               # noqa: E402
from s24 import stats_lib as ST                                          # noqa: E402

OUT = os.path.join(ROOT, "s26", "results", "i_identity_audit.json")
PINNED = ("peptide_folds.json", "peptide_clusters.json")
LANE_E_HASHES = os.path.join(ROOT, "s26", "results", "pinned_hashes.json")
#: the four self-copies the record names (state brief 3, ARCHITECTURE 2.1): (copy, carrier)
KNOWN_SELF_COPIES = (("1CEK", "1A11"), ("2FBU", "2LMF"), ("2P5H", "2P5J"), ("6B9K", "1U6V"))
N_FOLDS, SEED = 5, 0
THRESHOLD = D.IDENTITY_THRESHOLD


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def hashes():
    return {p: {"sha256": sha256(os.path.join(ROOT, p)),
                "bytes": os.path.getsize(os.path.join(ROOT, p))} for p in PINNED}


def containment_many(a, bs):
    """`identity` normalised by the SHORTER sequence, batched; the substring test is applied
    separately in `cluster` so its (redundant) contribution is counted, not assumed."""
    n = len(a)
    lens = np.array([len(b) for b in bs], float)
    return D._nw_matches(a, list(bs)) / np.maximum(np.minimum(lens, n), 1)


def cluster(seqs, mode):
    """Single-linkage clusters, sequence -> id by first appearance.

    `mode="longer"` reproduces `core.data.clusters` statement for statement (the admissible
    composition bound, then the exact alignment on the survivors).  `mode="shorter"` is the
    corrected measure: the bound and the alignment are both normalised by the shorter length,
    and a verbatim-substring pass follows.  The partition is unique whatever the union order,
    and ids are assigned by the lowest-index member, so the id map is reproducible.
    """
    C, lens = D.composition_matrix(seqs)
    parent = list(range(len(seqs)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    stats = {"n_bound_survivors": 0, "n_alignment_unions": 0, "n_substring_only_unions": 0}
    for a in range(len(seqs)):
        if mode == "longer":
            bound = np.asarray(D.max_possible_identity_many(seqs[a], C, lens), float).copy()
        else:
            shared = np.minimum(C[a][None, :], C).sum(1)
            bound = shared / np.maximum(np.minimum(lens, lens[a]), 1)
            bound = np.asarray(bound, float).copy()
        bound[:a + 1] = 0.0
        cand = [int(b) for b in np.flatnonzero(bound >= THRESHOLD) if find(a) != find(b)]
        if cand:
            stats["n_bound_survivors"] += len(cand)
            bs = [seqs[b] for b in cand]
            vals = (D.identity_many(seqs[a], bs) if mode == "longer"
                    else containment_many(seqs[a], bs))
            for b, v in zip(cand, vals):
                if v >= THRESHOLD and find(a) != find(b):
                    parent[find(b)] = find(a)
                    stats["n_alignment_unions"] += 1
        if mode == "shorter":
            sa = seqs[a]
            for b in range(a + 1, len(seqs)):
                if find(a) != find(b) and (sa in seqs[b] or seqs[b] in sa):
                    parent[find(b)] = find(a)
                    stats["n_substring_only_unions"] += 1
    roots, assign = {}, {}
    for i, s in enumerate(seqs):
        assign[s] = roots.setdefault(find(i), len(roots))
    return assign, stats


def fold_map(assign, n_folds=N_FOLDS, seed=SEED):
    """`core.data.folds`' shuffle, statement for statement, on an in-memory cluster map."""
    cl = sorted(set(assign.values()))
    rng = np.random.default_rng(seed)
    rng.shuffle(cl)
    fold_of = {c: i % n_folds for i, c in enumerate(cl)}
    return {s: fold_of[c] for s, c in assign.items()}


def members(assign, seq):
    cid = assign[seq]
    return [s for s, c in assign.items() if c == cid]


def main():
    t0 = time.time()
    before = hashes()
    lane_e = None
    if os.path.exists(LANE_E_HASHES):
        with open(LANE_E_HASHES) as fh:
            lane_e = json.load(fh)["files"]

    with open(os.path.join(ROOT, "peptide_clusters.json")) as fh:
        pinned_clusters = json.load(fh)
    with open(os.path.join(ROOT, "peptide_folds.json")) as fh:
        pinned_folds = json.load(fh)
    assert pinned_clusters.get("threshold") == THRESHOLD
    assert pinned_folds.get("n_folds") == N_FOLDS and pinned_folds.get("seed") == SEED
    pc, pf = pinned_clusters["assign"], pinned_folds["assign"]

    db = D.load()                                     # reads peptide_db.npz; never writes it
    seqs = [p.seq for p in db]
    pdb_of = {p.seq: p.pdb for p in db}
    print(f"database: {len(seqs)} peptides; pinned clusters {len(set(pc.values()))}, "
          f"pinned folds over {len(pf)} sequences")

    # ---- 1. the in-memory copy of the pinned convention must reproduce the pinned files
    t = time.time()
    a_long, st_long = cluster(seqs, "longer")
    f_long = fold_map(a_long)
    repro = {"clusters_identical": a_long == pc, "folds_identical": f_long == pf,
             "n_clusters": len(set(a_long.values())), "seconds": round(time.time() - t, 2),
             **st_long}
    print(f"reproduction of the pinned convention: clusters identical {repro['clusters_identical']}, "
          f"folds identical {repro['folds_identical']} ({repro['seconds']} s)")
    if not (repro["clusters_identical"] and repro["folds_identical"]):
        print("  !! the in-memory copy does NOT reproduce the pinned files; every number below "
              "is against the in-memory reproduction, not the pinned map")

    # ---- 2. the corrected measure
    t = time.time()
    a_short, st_short = cluster(seqs, "shorter")
    f_short = fold_map(a_short)
    corrected = {"n_clusters": len(set(a_short.values())), "seconds": round(time.time() - t, 2),
                 **st_short}
    changed_membership = [s for s in seqs if set(members(a_short, s)) != set(members(pc, s))]
    label_changed = [s for s in seqs if f_short[s] != pf[s]]
    corrected["n_sequences_whose_cluster_membership_changes"] = len(changed_membership)
    corrected["n_sequences_whose_naive_fold_label_changes"] = len(label_changed)
    print(f"corrected clustering: {corrected['n_clusters']} clusters (pinned "
          f"{repro['n_clusters']}); {len(changed_membership)} of {len(seqs)} sequences gain or "
          f"lose a cluster mate; a naive re-derivation of the folds would relabel "
          f"{len(label_changed)} of {len(seqs)} ({corrected['seconds']} s)")

    # ---- 3. the 126 dev targets (tuning126, public)
    from s12 import instrument as I
    T = I.targets()
    assert len(T) == 126
    dev = []
    for tg in T:
        s = tg["seq"]
        assert pf[s] == tg["fold"], (tg["pdb"], pf[s], tg["fold"])
        gained = sorted(set(members(a_short, s)) - set(members(pc, s)))
        lost = sorted(set(members(pc, s)) - set(members(a_short, s)))
        leak = [g for g in gained if pf[g] != pf[s]]      # a new mate sitting in another
        #                                                   pinned fold = inside this
        #                                                   target's fold model's training set
        row = {"pdb": tg["pdb"], "n": tg["n"], "pinned_fold": pf[s],
               "naive_fold_under_corrected": f_short[s],
               "label_changes": bool(f_short[s] != pf[s]),
               "n_new_mates": len(gained), "n_lost_mates": len(lost),
               "n_new_mates_in_other_pinned_folds": len(leak),
               "new_mate_lengths": [len(g) for g in gained],
               "new_mate_pinned_folds": [pf[g] for g in gained],
               "new_mate_containment": [round(float(D.identity(s, g, norm="shorter")), 4)
                                        for g in gained],
               "new_mate_longer_identity": [round(float(D.identity(s, g)), 4) for g in gained]}
        dev.append(row)
    dev_label = [r["pdb"] for r in dev if r["label_changes"]]
    dev_mates = [r["pdb"] for r in dev if r["n_new_mates"]]
    dev_leak = [r["pdb"] for r in dev if r["n_new_mates_in_other_pinned_folds"]]
    print(f"dev126: naive fold label changes on {len(dev_label)}/126; new cluster mates on "
          f"{len(dev_mates)}/126; new mates sitting in ANOTHER pinned fold (the leak) on "
          f"{len(dev_leak)}/126: {dev_leak}")

    # ---- 4. the four known self-copies
    known = []
    for copy, carrier in KNOWN_SELF_COPIES:
        pc_, pr_ = D.by_pdb(copy), D.by_pdb(carrier)
        row = {"copy": copy, "carrier": carrier,
               "copy_in_peptide_db": pc_ is not None, "carrier_in_peptide_db": pr_ is not None}
        if pc_ is not None and pr_ is not None:
            a, b = pc_.seq, pr_.seq
            row.update({
                "copy_len": len(a), "carrier_len": len(b),
                "verbatim_substring": bool(a in b or b in a),
                "longer_identity": round(float(D.identity(a, b)), 4),
                "shorter_identity": round(float(D.identity(a, b, norm="shorter")), 4),
                "pinned_same_cluster": pc[a] == pc[b],
                "pinned_folds": [pf[a], pf[b]],
                "corrected_same_cluster": a_short[a] == a_short[b],
                "copy_is_dev_target": copy in {t["pdb"] for t in T},
                "carrier_is_dev_target": carrier in {t["pdb"] for t in T},
            })
            row["would_be_moved_into_one_fold"] = bool(
                row["corrected_same_cluster"] and pf[a] != pf[b])
        else:
            row["note"] = ("the carrier is not a peptide-database entry, so it is outside the "
                           "identity clustering's reach (a fragment-bank window); the "
                           "corrected clustering cannot move it")
        known.append(row)
        print(f"  {copy} in {carrier}: {json.dumps(row)}")

    # ---- 5. the sealed benchmark: sequences only, counts only, through the permitted call
    bench = core.backend("data").benchmark()
    bseq = [p.seq for p in bench]
    b_label = sum(1 for s in bseq if f_short[s] != pf[s])
    b_leak = 0
    b_mates = 0
    for s in bseq:
        gained = set(members(a_short, s)) - set(members(pc, s))
        b_mates += bool(gained)
        b_leak += bool([g for g in gained if pf[g] != pf[s]])
    benchmark = {"n": len(bseq), "n_naive_fold_label_changes": int(b_label),
                 "n_with_new_cluster_mates": int(b_mates),
                 "n_with_new_mates_in_other_pinned_folds": int(b_leak),
                 "note": "counts only; names never read from the manifest by this script"}
    print(f"benchmark60 (counts only): naive label changes {b_label}/{len(bseq)}; new mates "
          f"{b_mates}/{len(bseq)}; new mates in another pinned fold {b_leak}/{len(bseq)}")

    # ---- 6. nothing was written to the pinned files
    after = hashes()
    untouched = before == after
    matches_lane_e = (lane_e is not None and all(
        lane_e.get(p, {}).get("sha256") == after[p]["sha256"] for p in PINNED))
    print(f"pinned files unchanged: {untouched}; match lane E's pinned_hashes.json: {matches_lane_e}")
    assert untouched, "the pinned files changed during the audit -- STOP"

    payload = {
        "label": "S26 lane I identity audit (in memory; pinned folds untouched)",
        "threshold": THRESHOLD, "n_folds": N_FOLDS, "seed": SEED,
        "n_sequences": len(seqs),
        "pinned_hashes_before": before, "pinned_hashes_after": after,
        "pinned_unchanged": untouched, "pinned_match_lane_E": matches_lane_e,
        "reproduction_of_pinned_convention": repro,
        "corrected": corrected,
        "dev126": {"n_naive_fold_label_changes": len(dev_label), "label_changes": dev_label,
                   "n_with_new_cluster_mates": len(dev_mates), "with_new_mates": dev_mates,
                   "n_with_new_mates_in_other_pinned_folds": len(dev_leak),
                   "with_new_mates_in_other_pinned_folds": dev_leak,
                   "rows": dev},
        "known_self_copies": known,
        "benchmark60": benchmark,
        "seconds_total": round(time.time() - t0, 1),
    }
    ST.save_atomic(OUT, payload, complete_keys=("pdb", "pinned_fold", "label_changes"),
                   rows=dev, n_expected=126, module_file=__file__)
    print(f"wrote {OUT} ({payload['seconds_total']} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
