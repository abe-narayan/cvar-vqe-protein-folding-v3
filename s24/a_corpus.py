"""SPRINT 24 / WORKSTREAM A / RUN 1 -- the corpus census and the exclusion audit.

Forks pre-registered in `s24/PREREG_A.md` and sent to the coordinator before this ran.
Summary of the six axes, so the operator choices travel with the artefact:

  functional     leak = {exact seq == , verbatim substring EITHER direction,
                 NW identity >= 0.6 normalised by the LONGER seq, source-PDB-ID equality}
                 NOT taken: containment >= 0.6 (measured to be AT the null for 9-16mers),
                 NOT taken: structural homology to benchmark natives (forbidden: needs
                 benchmark structures).
  basis          SOURCE CHAIN level, propagated to every window derived from it.
                 NOT taken: window-level only.
  readout        the exclusion list is the set of SOURCE IDS; reported with window counts
                 and per-target universe fractions.  NOT taken: a bare target count.
  normalisation  identity normalised by the LONGER sequence (the production rule).
                 NOT taken: normalising by the shorter (= containment).
  null           12 composition-matched random sequences per length, same criteria, same
                 bank, seeded from s15.seed.stable_rng.  NOT taken: no null.
  THE LABEL      "leak" = a corpus source chain that IS a held-out target, CONTAINS one
                 verbatim, is >=0.6 identical to one, or is a window of the SAME DEPOSIT.
                 NOT taken: "the generator's samples resemble the native" (circular).

WHAT THE CORPUS IS, established by reading the builder rather than by assumption
-------------------------------------------------------------------------------
`s8/generate.py:stage_univ` writes one npz per tuning target.  For target p in fold f:

    peps  = [q for q in peptide_db.load()  if folds[q.seq] != f and q.seq != p.seq]
    frags = distogram._fold_fragments(f, 5, threshold=0.6)
    W     = every length-n window of (peps + frags), IN THAT ORDER

so the window universe is a deterministic function of two banks -- `peptide_db` (787
chains, 8-26 aa, from pdbs/ + pdbs_ext/) and `fragment_db` (6,003 windows of 9-20 aa cut
from 1,001 deposits in prots/) -- and the per-window provenance is exactly reconstructible
without touching the 18k-row coordinate arrays.  This module reconstructs it and asserts
the reconstruction against every one of the 126 npz files.

Usage:  python -m s24.a_corpus            (census + audit + null + artefact)
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

import peptide_db as pdb                                                # noqa: E402
import fragment_db as fdb                                              # noqa: E402
import distogram as dgm                                                # noqa: E402
import core.data as cdata                                              # noqa: E402
from s7 import debias                                                  # noqa: E402
from s15.seed import stable_rng                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
UNIV = os.path.join(ROOT, "s8", "generate_univ")
os.makedirs(RESULTS, exist_ok=True)

IDENT_T = pdb.IDENTITY_THRESHOLD          # 0.6, the production rule
NULL_REPS = 12
NULL_LENGTHS = (9, 10, 11, 12, 13, 14, 15, 16)


# --------------------------------------------------------------------------- atomic io
def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, path)
    return path


def sha(obj):
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


# --------------------------------------------------------------------------- criteria
def strip_id(pid):
    """'1A6G_96' -> '1A6G';  '1ABC' -> '1ABC'.  Deposit id, chain/offset stripped."""
    return str(pid).split("_")[0].upper()[:4]


def verbatim(a, b):
    """True if either sequence occurs verbatim inside the other."""
    return (a in b) or (b in a)


def audit_pair(src_seq, tgt_seq, ct=None):
    """Every criterion for one (corpus chain, held-out target) pair."""
    out = {}
    out["exact"] = src_seq == tgt_seq
    out["substring"] = verbatim(src_seq, tgt_seq)
    ub = pdb.max_possible_identity(tgt_seq, src_seq, ct)
    out["identity"] = pdb.identity(tgt_seq, src_seq) if ub >= IDENT_T else 0.0
    out["homolog"] = out["identity"] >= IDENT_T
    return out


# --------------------------------------------------------------------------- census
def census():
    peps = list(pdb.load())
    frags = list(fdb.load())
    frag_src = sorted({strip_id(f.pdb) for f in frags})
    pep_src = sorted({strip_id(p.pdb) for p in peps})
    return peps, frags, pep_src, frag_src


def heldout_sets():
    tun = debias.tuning_targets()                 # the 126 dev instrument
    dev = cdata.dev_set(24)
    bench = cdata.benchmark()                     # IDENTIFIERS + SEQUENCES ONLY
    return {
        "tuning126": [(p.pdb, p.seq) for p in tun],
        "dev24": [(p.pdb, p.seq) for p in dev],
        "bench60": [(p.pdb, p.seq) for p in bench],
    }


# --------------------------------------------------------------- provenance of a universe
def universe_provenance(target, folds, frag_cache):
    """Exact (source id, kind, n_windows) list for one target's shipped universe.

    Reproduces `s8.generate._windows_all`'s iteration order: peptides first, then
    fragments, each contributing (len(seq) - n + 1) windows when len(seq) >= n.
    """
    f = folds[target.seq]
    peps = [q for q in pdb.load() if folds[q.seq] != f and q.seq != target.seq]
    frs = frag_cache[f]
    rows, tot = [], 0
    for q, kind in [(q, "peptide") for q in peps] + [(q, "fragment") for q in frs]:
        m = len(q.seq)
        if m < target.n:
            continue
        k = m - target.n + 1
        rows.append((q.pdb, q.seq, kind, tot, tot + k))
        tot += k
    return rows, tot


# --------------------------------------------------------------------------- main audit
def main():
    t0 = time.time()
    peps, frags, pep_src, frag_src = census()
    held = heldout_sets()
    folds = pdb.folds(5)

    print(f"corpus: {len(peps)} peptide chains ({len(pep_src)} deposits), "
          f"{len(frags)} fragments ({len(frag_src)} deposits)", flush=True)

    # ---- every corpus chain, with its identifier and sequence
    #
    # HAZARD, found here: `fragment_db` ids are `<PDB>_<start>` and are NOT UNIQUE -- the
    # same start yields one fragment per length in LENGTHS, all carrying the same id.  666
    # of the 6,003 ids collide.  Any downstream code that keys fragments by `f.pdb` silently
    # merges 9-mers with 20-mers.  The corpus primary key here is `kind:id:len`.
    corpus = ([dict(key=f"peptide:{p.pdb}:{p.n}", id=p.pdb, dep=strip_id(p.pdb),
                    seq=p.seq, kind="peptide") for p in peps]
              + [dict(key=f"fragment:{f.pdb}:{f.n}", id=f.pdb, dep=strip_id(f.pdb),
                      seq=f.seq, kind="fragment") for f in frags])
    assert len({c["key"] for c in corpus}) == len(corpus), "corpus key is not a key"
    pep_dep = {c["dep"] for c in corpus if c["kind"] == "peptide"}
    frg_dep = {c["dep"] for c in corpus if c["kind"] == "fragment"}
    print(f"deposit overlap peptide-bank vs fragment-bank: {len(pep_dep & frg_dep)}"
          f"  ({sorted(pep_dep & frg_dep)[:10]})", flush=True)

    # ---- audit every corpus chain against every held-out target
    tgt_comp = {}
    for setname, rows in held.items():
        for pid, seq in rows:
            tgt_comp[(setname, pid)] = (seq, pdb._composition(seq))

    hits = {}          # corpus id -> list of hit records
    for c in corpus:
        for (setname, pid), (tseq, ct) in tgt_comp.items():
            r = audit_pair(c["seq"], tseq, ct)
            dep = c["dep"] == strip_id(pid)
            if not (r["exact"] or r["substring"] or r["homolog"] or dep):
                continue
            hits.setdefault(c["key"], []).append(dict(
                set=setname, target=pid, kind=c["kind"], dep=c["dep"],
                deposit_match=bool(dep), **r))
    print(f"audit done ({time.time()-t0:.0f}s): {len(hits)} corpus chains hit", flush=True)

    # ---- per-criterion counts (published even when zero)
    crits = ["exact", "substring", "homolog", "deposit_match"]
    by_crit = {s: {c: set() for c in crits} for s in held}
    tgts_hit = {s: {c: set() for c in crits} for s in held}
    for cid, rs in hits.items():
        for r in rs:
            for c in crits:
                if r[c]:
                    by_crit[r["set"]][c].add(cid)
                    tgts_hit[r["set"]][c].add(r["target"])
    crit_table = {s: {c: dict(n_corpus_chains=len(by_crit[s][c]),
                              n_targets_hit=len(tgts_hit[s][c]),
                              targets=sorted(tgts_hit[s][c]))
                      for c in crits} for s in held}

    # ---- kind split, for H-A2 (peptide half vs fragment half)
    kindsplit = {}
    for s in held:
        d = {}
        for k in ("peptide", "fragment"):
            n_tot = len(peps) if k == "peptide" else len(frags)
            ids = {cid for cid, rs in hits.items()
                   if any(r["set"] == s and r["kind"] == k for r in rs)}
            d[k] = dict(n_hit=len(ids), n_total=n_tot, rate=len(ids) / n_tot)
        kindsplit[s] = d

    # ---- THE EXCLUSION LIST.
    # Pre-registered BASIS fork: the verdict is taken at SOURCE-CHAIN level and PROPAGATED
    # TO EVERY WINDOW DERIVED FROM IT.  For a peptide the source chain is the entry; for a
    # fragment the source chain is its PARENT DEPOSIT, so one contaminated fragment removes
    # every fragment cut from that deposit.  That is the strict reading and it is the one
    # shipped.  The looser chain-only count is reported beside it so the operator is visible.
    hit_chains = sorted(hits)
    hit_deps = {r["dep"] for rs in hits.values() for r in rs}
    exclude = sorted(c["key"] for c in corpus if c["dep"] in hit_deps)
    exclude_dep = sorted(hit_deps)
    permitted = sorted(c["key"] for c in corpus if c["dep"] not in hit_deps)

    # ---- per-target universe impact, and the reconstruction assertion
    frag_cache = {f: list(dgm._fold_fragments(f, 5, IDENT_T)) for f in range(5)}
    print("fragment sets per fold: " +
          ", ".join(f"{f}:{len(frag_cache[f])}" for f in range(5)), flush=True)
    tun = {p.pdb: p for p in debias.tuning_targets()}
    ex = set(hit_deps)
    per_target, recon_ok = {}, True
    for pid in sorted(tun):
        rows, tot = universe_provenance(tun[pid], folds, frag_cache)
        z = np.load(os.path.join(UNIV, f"{pid}.npz"), allow_pickle=True)
        nw = int(z["W"].shape[0])
        ok = (nw == tot)
        recon_ok &= ok
        bad = sum(b - a for sid, sq, kd, a, b in rows if strip_id(sid) in ex)
        badp = sum(b - a for sid, sq, kd, a, b in rows
                   if strip_id(sid) in ex and kd == "peptide")
        per_target[pid] = dict(n_windows=nw, n_windows_reconstructed=tot,
                               reconstruction_ok=bool(ok),
                               n_sources=len(rows), n_excluded_windows=int(bad),
                               n_excluded_windows_peptide=int(badp),
                               frac_excluded=float(bad) / max(nw, 1))
    fr = np.array([v["frac_excluded"] for v in per_target.values()])
    print(f"reconstruction exact on all 126: {recon_ok}; "
          f"excluded window fraction mean {fr.mean():.4f} max {fr.max():.4f}", flush=True)

    # ---- the null: composition-matched random sequences, same criteria, same bank
    null = {}
    aa = np.array(list("ARNDCQEGHILKMFPSTWYV"))
    allseq = "".join(c["seq"] for c in corpus)
    p_aa = np.array([allseq.count(a) for a in aa], float)
    p_aa /= p_aa.sum()
    for n in NULL_LENGTHS:
        rng = stable_rng("s24A", "null", n)
        ident, sub, exa = [], 0, 0
        for _ in range(NULL_REPS):
            s = "".join(rng.choice(aa, size=n, p=p_aa))
            cs = pdb._composition(s)
            best = 0.0
            for c in corpus:
                if pdb.max_possible_identity(s, c["seq"], cs) <= best:
                    continue
                best = max(best, pdb.identity(s, c["seq"]))
            ident.append(best)
            sub += any(verbatim(s, c["seq"]) for c in corpus)
            exa += any(s == c["seq"] for c in corpus)
        null[n] = dict(max_identity_mean=float(np.mean(ident)),
                       max_identity_max=float(np.max(ident)),
                       n_with_verbatim_substring=int(sub), n_exact=int(exa),
                       n_reps=NULL_REPS)
        print(f"null n={n}: max-identity mean {np.mean(ident):.3f} "
              f"max {np.max(ident):.3f}, substring hits {sub}/{NULL_REPS}", flush=True)

    # ---- artefact
    art = dict(
        corpus_definition=dict(
            builder="s8/generate.py:stage_univ -> _windows_all(peps + frags)",
            peptide_bank=dict(module="peptide_db", cache="peptide_db.npz",
                              n_chains=len(peps), n_deposits=len(pep_src),
                              len_range=[8, 26], sources=["pdbs/", "pdbs_ext/"]),
            fragment_bank=dict(module="fragment_db", cache="fragment_db.npz",
                               n_fragments=len(frags), n_deposits=len(frag_src),
                               len_range=[9, 20], stride=5, source="prots/"),
            fold_rule="peptides: folds[seq] != target fold AND seq != target seq; "
                      "fragments: distogram._fold_fragments(fold, 5, 0.6)",
        ),
        criteria=dict(identity_threshold=IDENT_T,
                      identity_normalisation="longer sequence (peptide_db.identity)",
                      substring="verbatim, either direction",
                      deposit="4-character PDB id equality after stripping _offset"),
        heldout=dict((k, [p for p, _ in v]) for k, v in held.items()),
        criterion_table=crit_table,
        kind_split=kindsplit,
        null=null,
        n_corpus_chains=len(corpus),
        n_hit_chains=len(hit_chains),
        n_excluded_chains=len(exclude),
        n_excluded_deposits=len(exclude_dep),
        n_permitted_chains=len(permitted),
        deposit_overlap_peptide_fragment=sorted(pep_dep & frg_dep),
        fragment_id_collisions=len(frags) - len({f["id"] for f in corpus
                                                 if f["kind"] == "fragment"}),
        reconstruction_exact_on_126=bool(recon_ok),
        per_target=per_target,
        hits=hits,
    )
    art["exclusion_list"] = exclude
    art["excluded_deposits"] = exclude_dep
    art["permitted_chains"] = permitted
    key = dict(exclusion_list=exclude, permitted_chains=permitted,
               criteria=art["criteria"], heldout=art["heldout"])
    art["CORPUS_HASH"] = sha(key)[:16]
    art["complete"] = bool(
        recon_ok
        and set(per_target) == set(tun)
        and all(set(crit_table[s]) == set(crits) for s in held)
        and set(null) == set(NULL_LENGTHS)
        and len(exclude) + len(permitted) == len(corpus))

    write_json(os.path.join(RESULTS, "a_corpus_audit.json"), art)
    write_json(os.path.join(RESULTS, f"a_corpus_permitted_{art['CORPUS_HASH']}.json"),
               dict(CORPUS_HASH=art["CORPUS_HASH"], criteria=art["criteria"],
                    heldout=art["heldout"], exclusion_list=exclude,
                    excluded_deposits=exclude_dep, permitted_chains=permitted,
                    n_permitted=len(permitted), n_excluded=len(exclude),
                    complete=art["complete"]))
    print(json.dumps({k: art[k] for k in
                      ("n_corpus_chains", "n_excluded_chains", "n_excluded_deposits",
                       "n_permitted_chains", "reconstruction_exact_on_126",
                       "CORPUS_HASH", "complete")}, indent=1))
    print(json.dumps(crit_table, indent=1))
    print(json.dumps(kindsplit, indent=1))
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
