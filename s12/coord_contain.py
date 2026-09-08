"""COORDINATOR AUDIT 1 -- the TRAINING-SET containment leak, which has never been measured.

The record prices a RETRIEVAL leak (a >=0.6-identity window inside the K=500 pool) at
+0.0004 A on tuning126.  That is not the same object as the leak hazard H3 describes.

`peptide_db.identity` normalises the Needleman-Wunsch match count by the LONGER sequence.
`clusters()` thresholds it at 0.6 and `folds()` partitions those clusters.  So a 13-residue
target sitting VERBATIM inside a 25-residue library member scores 13/25 = 0.52, lands in a
DIFFERENT identity cluster, and can therefore land in a DIFFERENT FOLD -- which means that
member is in the TRAINING SET of the distogram that scores the target.  The distogram then
learned the target's own pair distances from a chain that contains it.

This module measures that pathway directly, for the 126 tuning targets, over BOTH training
corpora the distogram uses (787 peptides + the fold's fragments), using `containment`
(normalised by the SHORTER sequence) -- the metric H3 says is the correct one.

Outputs s12/results/coord_containment.json.  Nothing here changes any pinned object.
"""
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from s12 import instrument as I           # noqa: E402
import peptide_db as pdb                  # noqa: E402
import fragment_db as fdb                 # noqa: E402


def containment(a: str, b: str, gap: float = -1.0) -> float:
    """NW match count normalised by the SHORTER sequence (core.data.containment)."""
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return 0.0
    prev = [(gap * j, 0) for j in range(m + 1)]
    for i in range(1, n + 1):
        cur = [(gap * i, 0)]
        ai = a[i - 1]
        for j in range(1, m + 1):
            match = 1.0 if ai == b[j - 1] else 0.0
            diag = (prev[j - 1][0] + match, prev[j - 1][1] + int(match > 0))
            up = (prev[j][0] + gap, prev[j][1])
            left = (cur[j - 1][0] + gap, cur[j - 1][1])
            best = diag
            if up[0] > best[0]:
                best = up
            if left[0] > best[0]:
                best = left
            cur.append(best)
        prev = cur
    return prev[m][1] / max(min(n, m), 1)


def substring_hit(t: str, s: str) -> bool:
    return t in s


def main():
    tg = I.targets()
    folds = pdb.folds(5)
    peps = list(pdb.load())
    frags = list(fdb.load(False))
    print(f"{len(tg)} targets, {len(peps)} peptides, {len(frags)} fragments", flush=True)

    # cheap pre-screen: 3-mer overlap, then exact substring, then full containment
    rows = []
    for t in tg:
        seq, f = t["seq"], t["fold"]
        tk = {seq[i:i + 3] for i in range(len(seq) - 2)}
        hits = {"peptide_train": [], "peptide_heldout": [], "fragment_train": []}
        for q in peps:
            if q.seq == seq:
                continue
            if not (tk & {q.seq[i:i + 3] for i in range(len(q.seq) - 2)}):
                continue
            c = containment(seq, q.seq)
            if c < 0.6:
                continue
            ident = pdb.identity(seq, q.seq)
            in_train = folds[q.seq] != f          # trained the distogram that scores `seq`
            rec = {"pdb": q.pdb, "seq": q.seq, "containment": round(c, 4),
                   "identity": round(ident, 4), "verbatim": substring_hit(seq, q.seq),
                   "fold": int(folds[q.seq]), "in_distogram_training_set": bool(in_train)}
            hits["peptide_train" if in_train else "peptide_heldout"].append(rec)
        # fragments: only this fold's fragment set trains the model, and `_fold_fragments`
        # already drops fragments above 0.6 IDENTITY to any peptide of the fold.  The
        # containment metric is the one that can slip through.
        for q in frags:
            if not (tk & {q.seq[i:i + 3] for i in range(len(q.seq) - 2)}):
                continue
            c = containment(seq, q.seq)
            if c < 0.6:
                continue
            hits["fragment_train"].append(
                {"pdb": q.pdb, "seq": q.seq, "containment": round(c, 4),
                 "identity": round(pdb.identity(seq, q.seq), 4),
                 "verbatim": substring_hit(seq, q.seq)})
        rows.append({"pdb": t["pdb"], "seq": seq, "n": t["n"], "fold": f,
                     "n_peptide_train": len(hits["peptide_train"]),
                     "n_peptide_heldout": len(hits["peptide_heldout"]),
                     "n_fragment_train": len(hits["fragment_train"]),
                     "verbatim_in_training": any(h["verbatim"] for h in hits["peptide_train"])
                                             or any(h["verbatim"] for h in hits["fragment_train"]),
                     "hits": hits})
        print(f"  {t['pdb']:5s} n={t['n']:2d} f={f}  train_pep={len(hits['peptide_train']):2d} "
              f"held={len(hits['peptide_heldout']):2d} frag={len(hits['fragment_train']):3d} "
              f"verbatim={rows[-1]['verbatim_in_training']}", flush=True)

    n_leak = sum(r["n_peptide_train"] > 0 or r["n_fragment_train"] > 0 for r in rows)
    n_verb = sum(r["verbatim_in_training"] for r in rows)
    out = {"what": "TRAINING-SET containment leak into the leave-fold-out distogram, 126 tuning targets",
           "metric": "Needleman-Wunsch match count / len(shorter)  (core.data.containment); "
                     "the production fold split uses /len(longer), which is the hazard",
           "n_targets": len(rows),
           "n_targets_with_any_training_containment_ge_0.6": n_leak,
           "n_targets_with_verbatim_copy_in_training_set": n_verb,
           "note": "in_distogram_training_set is True when the library member's fold differs "
                   "from the target's fold, i.e. it was used to fit the model that scores it",
           "per_target": rows}
    p = I.write("coord_containment", out)
    print(f"\n{n_leak}/126 targets have a >=0.6-containment member in their distogram's "
          f"training set; {n_verb}/126 have a VERBATIM copy.\nwrote {p}")


if __name__ == "__main__":
    main()
