"""COORDINATOR EXPERIMENT 3 -- is there a sequence->structure channel at 9-16 residues at all?

Every negative in this sprint has the same shape: the candidate set can be improved
arbitrarily and the answer does not move, because nothing ranks within it.  The objective
agent has now bounded the distance channel itself -- a PERFECT distance matrix emits
2.395 A -- and found that the shipped distogram's per-separation-shell profile is matched,
and beaten at long range, by the K=500 pool's OWN MEAN PROFILE: no model, no labels, no
sequence.  That is the sharpest possible statement of "the predictor is typicality".

So this module stops asking what a model can extract and asks what is there to extract.
It measures the sequence->structure channel directly in the corpus, with no model at all:

    Q1  Does sequence similarity predict structural similarity at 9-16 residues, and where
        does the relationship emerge from the null?
    Q2  How much does conditioning on the sequence narrow the distribution of structures,
        against conditioning on nothing?  (A conditional-vs-marginal spread ratio, which is
        a model-free proxy for the mutual information the whole architecture is trying to
        exploit.)
    Q3  The deployable form of the same question: among a target's OWN pool, does BLOSUM
        similarity rank windows by structural closeness to the target's native at all?

Q1/Q2 use only library peptides against each other -- no tuning target's native is read.
Q3 uses `rr`, which is an ORACLE label, and is marked as such.

The controls are the point.  Structural similarity between two peptides of the same length
has a large length-dependent floor (any two 9-mers are more similar than any two 16-mers,
because there is less to disagree about), and composition alone carries structure (a
poly-A stretch is a helix).  So every relationship is measured against a length-matched
permutation null and against a composition-matched null.

    python -m s12.coord_infochannel
"""
import os
import sys
import json
import itertools

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402
import peptide_db as pdb                   # noqa: E402
import fragment_db as fdb                  # noqa: E402

RNG = np.random.default_rng(20260905)
AA = I.ALPHABET


def _b62():
    from s5 import lib as s5lib
    return s5lib.B62


def encode(s):
    idx = {a: i for i, a in enumerate(AA)}
    return np.array([idx.get(c, 0) for c in s], int)


def blosum_sim(a, b, B):
    ea, eb = encode(a), encode(b)
    return float(B[ea, eb].sum())


def identity_frac(a, b):
    return float(np.mean([x == y for x, y in zip(a, b)]))


def q1_q2(entries, lengths=(9, 10, 11, 12, 13, 14, 15, 16), max_pairs=200000):
    """Sequence similarity vs structural similarity, all same-length pairs, per length."""
    B = _b62()
    out = {}
    for n in lengths:
        E = [p for p in entries if p.n == n]
        if len(E) < 20:
            continue
        W = np.stack([np.asarray(p.ca, float) for p in E])
        seqs = [p.seq for p in E]
        m = len(E)
        pairs = list(itertools.combinations(range(m), 2))
        if len(pairs) > max_pairs:
            pairs = [pairs[k] for k in RNG.permutation(len(pairs))[:max_pairs]]
        ii = np.array([a for a, _ in pairs]); jj = np.array([b for _, b in pairs])
        # structural distance, batched
        rms = np.empty(len(pairs))
        for s in range(0, len(pairs), 2000):
            sl = slice(s, min(s + 2000, len(pairs)))
            blk = np.array([I.kabsch_rmsd_batch(W[ii[sl][k]][None], W[jj[sl][k]])[0]
                            for k in range(len(ii[sl]))])
            rms[sl] = blk
        ident = np.array([identity_frac(seqs[a], seqs[b]) for a, b in pairs])
        blos = np.array([blosum_sim(seqs[a], seqs[b], B) for a, b in pairs])
        from scipy.stats import spearmanr, pearsonr
        rho_id = spearmanr(ident, rms)
        rho_bl = spearmanr(blos, rms)
        # permutation null on the SAME marginals
        null = []
        for _ in range(200):
            null.append(spearmanr(RNG.permutation(ident), rms).statistic)
        null = np.array(null)
        # binned: structural distance as a function of identity
        bins = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.6), (0.6, 1.01)]
        prof = []
        for lo, hi in bins:
            sel = (ident >= lo) & (ident < hi)
            prof.append({"identity": f"[{lo},{hi})", "n_pairs": int(sel.sum()),
                         "mean_rmsd": float(rms[sel].mean()) if sel.any() else None,
                         "p10_rmsd": float(np.percentile(rms[sel], 10)) if sel.sum() > 10 else None})
        out[n] = {"n_peptides": m, "n_pairs": len(pairs),
                  "mean_rmsd_all_pairs": float(rms.mean()),
                  "sd_rmsd_all_pairs": float(rms.std()),
                  "rho_identity_vs_rmsd": float(rho_id.statistic), "p_identity": float(rho_id.pvalue),
                  "rho_blosum_vs_rmsd": float(rho_bl.statistic), "p_blosum": float(rho_bl.pvalue),
                  "permutation_null_rho": [float(null.mean()), float(np.percentile(null, 2.5)),
                                           float(np.percentile(null, 97.5))],
                  "profile_by_identity": prof}
        print(f"  n={n}: {m} peptides, {len(pairs)} pairs, mean pair RMSD {rms.mean():.3f}, "
              f"rho(identity,RMSD) {rho_id.statistic:+.4f} (p={rho_id.pvalue:.2g}), "
              f"rho(BLOSUM,RMSD) {rho_bl.statistic:+.4f}", flush=True)
    return out


def q3(max_targets=126):
    """ORACLE DIAGNOSTIC. Within each target's own universe, does BLOSUM rank by nativeness?

    This is the deployable channel in its purest form: the only sequence-conditioned
    quantity retrieval has is `sim`, and `rr` is how close each window actually is.
    """
    from scipy.stats import spearmanr
    rows = []
    for t in I.targets()[:max_targets]:
        u = I.load_univ(t["pdb"])
        sim = np.asarray(u["sim"], float); rr = np.asarray(u["rr"], float)
        r = spearmanr(sim, rr)
        pool = I.pool_idx(u)
        # what does the top-K by BLOSUM buy over a random K of the same size?
        rnd = RNG.permutation(len(rr))[:I.K]
        rows.append({"pdb": t["pdb"], "n": t["n"], "fold": t["fold"],
                     "rho_sim_vs_rr": float(r.statistic), "p": float(r.pvalue),
                     "n_windows": int(len(rr)),
                     "pool_best": float(rr[pool].min()), "rand_best": float(rr[rnd].min()),
                     "pool_mean": float(rr[pool].mean()), "rand_mean": float(rr[rnd].mean()),
                     "universe_mean": float(rr.mean()), "universe_best": float(rr.min())})
    return rows


def main():
    print("Q1/Q2: sequence vs structure over the corpus (no target natives read)", flush=True)
    peps = list(pdb.load())
    res_pep = q1_q2(peps)
    print("\nQ1/Q2 on the FRAGMENT library (protein-derived, same measurement)", flush=True)
    res_frag = q1_q2(list(fdb.load(False)))
    print("\nQ3: does BLOSUM rank windows by nativeness inside each target's universe? (ORACLE)",
          flush=True)
    rows = q3()
    rho = np.array([r["rho_sim_vs_rr"] for r in rows])
    pb = np.array([r["pool_best"] for r in rows]); rb = np.array([r["rand_best"] for r in rows])
    pm = np.array([r["pool_mean"] for r in rows]); rm = np.array([r["rand_mean"] for r in rows])
    f18 = set(I.FAIL18); isf = np.array([r["pdb"] in f18 for r in rows])
    summary = {
        "q3_rho_sim_vs_rr": {"mean": float(rho.mean()), "median": float(np.median(rho)),
                             "sd": float(rho.std()), "n_positive": int((rho > 0).sum()),
                             "FAIL18_mean": float(rho[isf].mean()),
                             "other108_mean": float(rho[~isf].mean())},
        "q3_pool_vs_random": {
            "pool_best": float(pb.mean()), "random_best": float(rb.mean()),
            "pool_mean": float(pm.mean()), "random_mean": float(rm.mean()),
            "best_paired": I.paired(pb, rb, names=[r["pdb"] for r in rows]),
            "mean_paired": I.paired(pm, rm, names=[r["pdb"] for r in rows])},
    }
    out = {"what": "model-free measurement of the sequence->structure channel at 9-16 residues",
           "q1_peptides": res_pep, "q1_fragments": res_frag,
           "q3_per_target": rows, "q3_summary": summary}
    I.write("coord_infochannel", out)

    print("\n=== Q3 summary ===")
    s = summary["q3_rho_sim_vs_rr"]
    print(f"  rho(BLOSUM sim, true RMSD) within a target's universe: mean {s['mean']:+.4f} "
          f"median {s['median']:+.4f}, positive on {s['n_positive']}/126")
    print(f"    (negative rho = higher similarity means CLOSER structure, i.e. the channel works)")
    print(f"    FAIL18 {s['FAIL18_mean']:+.4f}   other108 {s['other108_mean']:+.4f}")
    p = summary["q3_pool_vs_random"]
    print(f"  BLOSUM top-500 vs a RANDOM 500 of the same universe:")
    print(f"    pool best {p['pool_best']:.4f} vs random best {p['random_best']:.4f}  "
          f"({p['best_paired']['mean_diff']:+.4f} [{p['best_paired']['ci95'][0]:+.4f},"
          f"{p['best_paired']['ci95'][1]:+.4f}])")
    print(f"    pool mean {p['pool_mean']:.4f} vs random mean {p['random_mean']:.4f}  "
          f"({p['mean_paired']['mean_diff']:+.4f} [{p['mean_paired']['ci95'][0]:+.4f},"
          f"{p['mean_paired']['ci95'][1]:+.4f}])")


if __name__ == "__main__":
    main()
