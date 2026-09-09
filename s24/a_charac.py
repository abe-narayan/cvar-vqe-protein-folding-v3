"""SPRINT 24 / WORKSTREAM A / RUN 2 -- corpus characterisation and Lane C's decisive question.

Forks pre-registered in `s24/PREREG_A.md`, sent to and approved by the coordinator before
this ran, and reproduced here so the operator choices travel with the artefact:

  functional     near-duplicate collapse under Kabsch CA-RMSD, whole curve tau in
                 {0.5,1.0,1.5,2.0} A.  NOT taken as primary: torsion-space distance (the
                 endpoint is CA-RMSD and the control must match the operator's space);
                 reported as a SECONDARY and disagreement is reported.
  basis          window level, LENGTH-STRATIFIED 9..16, never pooled across lengths.
                 NOT taken: pooling lengths, which invents similarity from length alone.
  readout        effective independent count = greedy leader count at each tau, plus the
                 exp(H) effective sample size of the cluster-size distribution.
                 NOT taken: the raw window count, the number the corpus advertises.
  normalisation  per length, on the PERMITTED corpus only, reported absolutely AND as a
                 ratio to a plausible generator's parameter count.
                 NOT taken: an absolute count with no model-size referent.
  null           two PLAUSIBLE zero-information controls in the operator's space:
                 (i) a constant alpha-helix, (ii) torsions drawn from the per-residue-type
                 marginal Ramachandran table `s8/generate_rama.npz`.
                 NOT taken: uniform-on-the-torus (a WORSE measure, not an uninformative one).
  THE LABEL      "the corpus supports a from-scratch generator" = >=10x as many tau=1.0 A
                 independent windows per length band as the smallest plausible generator's
                 parameter count, AND genuine per-residue-type phi/psi multimodality rather
                 than one dominant basin.  Both halves stated in advance.
                 NOT taken: "supports a generator" = "there is a lot of data".

Pre-registered NO conditions, stated before the run so a NO is not a post-hoc rescue:
fewer than ~1e4 independent windows in ANY length band, OR >80% of windows in one
Ramachandran basin class per residue type, OR an independent count that does not clearly
exceed the plausible zero-information control.

The greedy-leader count is measured on SUBSAMPLES with a GROWTH CURVE
(N = 500,1000,2000,5000,10000), because the decisive question is not the count at one N --
it is whether the count SATURATES.  A saturating curve is "a small number of folds
repeated"; a curve still growing near-linearly at the largest N is genuinely new structure
per window.  This is stated here because it is the reading, not a convenience.

Usage:  python -m s24.a_charac cheap     # census, torsions, SS, vocabulary; no lock
        python -m s24.a_charac heavy     # the RMSD leader curves; TAKES LOCK_TRAIN
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import peptide_db as pdb                                                # noqa: E402
import fragment_db as fdb                                              # noqa: E402
from s15.seed import stable_rng                                        # noqa: E402
from s12 import instrument as I                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
LOCK = os.path.join(RESULTS, "LOCK_TRAIN")
os.makedirs(RESULTS, exist_ok=True)

HASH = "29e3b67e8ca0c03d"
PERMITTED = os.path.join(RESULTS, f"a_corpus_permitted_{HASH}.json")
LENGTHS = (9, 10, 11, 12, 13, 14, 15, 16)
TAUS = (0.5, 1.0, 1.5, 2.0)
NGRID = (500, 1000, 2000, 5000, 10000)
CTRL_NGRID = (500, 1000, 2000, 5000)      # controls capped; compared at MATCHED N only
ALPHABET = "ARNDCQEGHILKMFPSTWYV"


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, path)


# --------------------------------------------------------------------------- corpus
def permitted_chains():
    """The permitted source chains, as (key, seq, ca, phi, psi, kind, deposit)."""
    with open(PERMITTED) as fh:
        art = json.load(fh)
    keep = set(art["permitted_chains"])
    out = []
    for p in pdb.load():
        k = f"peptide:{p.pdb}:{p.n}"
        if k in keep:
            out.append((k, p.seq, np.asarray(p.ca, float), np.asarray(p.phi, float),
                        np.asarray(p.psi, float), "peptide", p.pdb.split("_")[0][:4]))
    for f in fdb.load():
        k = f"fragment:{f.pdb}:{f.n}"
        if k in keep:
            out.append((k, f.seq, np.asarray(f.ca, float), np.asarray(f.phi, float),
                        np.asarray(f.psi, float), "fragment", f.pdb.split("_")[0][:4]))
    assert len(out) == len(keep), f"{len(out)} != {len(keep)}"
    return out, art["CORPUS_HASH"]


def windows_of(chains, L):
    """Every length-L window of the permitted corpus.  Returns CA, PHI, PSI, SEQ, SRC."""
    ca, ph, ps, sq, src = [], [], [], [], []
    for idx, (key, seq, c, f, y, kind, dep) in enumerate(chains):
        m = len(seq)
        if m < L:
            continue
        for s in range(m - L + 1):
            ca.append(c[s:s + L]); ph.append(f[s:s + L]); ps.append(y[s:s + L])
            sq.append(seq[s:s + L]); src.append(idx)
    if not ca:
        return None
    return (np.stack(ca), np.stack(ph), np.stack(ps), sq, np.array(src, int))


# --------------------------------------------------------------------------- basins
def basin(phi, psi):
    """Coarse Ramachandran class.  0 alphaR, 1 beta/PPII, 2 alphaL, 3 other."""
    d = np.degrees
    p, q = d(phi), d(psi)
    out = np.full(p.shape, 3, np.int8)
    out[(p < 0) & (q > -120) & (q < 50)] = 0
    out[(p < 0) & ((q >= 50) | (q <= -120))] = 1
    out[(p >= 0) & (q > -60) & (q < 120)] = 2
    return out


# --------------------------------------------------------------------------- cheap
def cheap():
    t0 = time.time()
    chains, h = permitted_chains()
    assert h == HASH
    deps = sorted({c[6] for c in chains})
    kinds = {k: sum(1 for c in chains if c[5] == k) for k in ("peptide", "fragment")}
    print(f"permitted: {len(chains)} chains, {len(deps)} deposits, {kinds}", flush=True)

    out = dict(CORPUS_HASH=h, n_chains=len(chains), n_deposits=len(deps), kinds=kinds)

    # ---- length distribution of the SOURCE chains and of the WINDOWS
    out["chain_len_hist"] = {int(L): int(sum(1 for c in chains if len(c[1]) == L))
                             for L in sorted({len(c[1]) for c in chains})}
    wcount, wsrc, wdep = {}, {}, {}
    for L in LENGTHS:
        w = windows_of(chains, L)
        wcount[L] = 0 if w is None else int(len(w[0]))
        wsrc[L] = 0 if w is None else int(len(set(w[4].tolist())))
        wdep[L] = 0 if w is None else int(len({chains[i][6] for i in set(w[4].tolist())}))
    out["windows_per_length"] = wcount
    out["source_chains_per_length"] = wsrc
    out["source_deposits_per_length"] = wdep
    print("windows/length " + json.dumps(wcount), flush=True)

    # ---- residue vocabulary over the whole permitted corpus
    allseq = "".join(c[1] for c in chains)
    comp = {a: allseq.count(a) for a in ALPHABET}
    tot = sum(comp.values())
    out["residue_vocabulary"] = {a: comp[a] / tot for a in ALPHABET}
    out["n_residues"] = tot
    out["vocab_entropy_bits"] = float(-sum((v / tot) * math.log2(v / tot)
                                           for v in comp.values() if v))
    out["unknown_residue_chars"] = sorted(set(allseq) - set(ALPHABET))

    # ---- torsion marginals and per-residue-type basin composition (THE LABEL, half 2)
    PH = np.concatenate([c[3] for c in chains])
    PS = np.concatenate([c[4] for c in chains])
    AA = np.concatenate([np.array([ALPHABET.index(x) if x in ALPHABET else -1
                                   for x in c[1]]) for c in chains])
    ok = np.isfinite(PH) & np.isfinite(PS) & (AA >= 0)
    PH, PS, AA = PH[ok], PS[ok], AA[ok]
    B = basin(PH, PS)
    glob = np.bincount(B, minlength=4) / len(B)
    out["basin_global"] = dict(alphaR=float(glob[0]), beta=float(glob[1]),
                               alphaL=float(glob[2]), other=float(glob[3]))
    per = {}
    for a, ch in enumerate(ALPHABET):
        m = AA == a
        if m.sum() < 50:
            continue
        f = np.bincount(B[m], minlength=4) / m.sum()
        per[ch] = dict(n=int(m.sum()), alphaR=float(f[0]), beta=float(f[1]),
                       alphaL=float(f[2]), other=float(f[3]), max_frac=float(f.max()),
                       phi_mean_deg=float(np.degrees(np.arctan2(
                           np.sin(PH[m]).mean(), np.cos(PH[m]).mean()))),
                       psi_circ_sd_deg=float(np.degrees(math.sqrt(-2 * math.log(
                           max(1e-12, math.hypot(np.sin(PS[m]).mean(),
                                                 np.cos(PS[m]).mean())))))))
    out["basin_per_residue"] = per
    out["n_residues_with_max_basin_over_80pct"] = int(
        sum(1 for v in per.values() if v["max_frac"] > 0.80))
    out["n_residue_types_measured"] = len(per)
    print(f"global basins {out['basin_global']}; residue types with a >80% dominant "
          f"basin: {out['n_residues_with_max_basin_over_80pct']}/{len(per)}", flush=True)

    # ---- secondary-structure composition, per SOURCE CHAIN (exact, then sliced)
    from core import geometry as geo
    ssc = {"H": 0, "E": 0, "C": 0}
    ss_by_chain = {}
    bad = 0
    for key, seq, c, f, y, kind, dep in chains:
        try:
            bb = geo.build_backbone(f, y)
            s = geo.assign_secondary_structure(bb)
        except Exception:
            bad += 1
            continue
        ss_by_chain[key] = s
        for ch in s:
            if ch in ssc:
                ssc[ch] += 1
    n = sum(ssc.values())
    out["ss_composition"] = {k: v / max(n, 1) for k, v in ssc.items()}
    out["ss_failed_chains"] = bad
    print(f"SS composition {out['ss_composition']} (failed {bad})", flush=True)

    # ---- window-level SS-string diversity, per length: how many DISTINCT H/E/C strings
    ssdiv = {}
    for L in LENGTHS:
        strs = {}
        for key, seq, c, f, y, kind, dep in chains:
            s = ss_by_chain.get(key)
            if s is None or len(s) < L:
                continue
            for a in range(len(s) - L + 1):
                strs[s[a:a + L]] = strs.get(s[a:a + L], 0) + 1
        tot = sum(strs.values())
        p = np.array(list(strs.values()), float) / max(tot, 1)
        ssdiv[L] = dict(n_windows=int(tot), n_distinct=len(strs),
                        n_possible=int(3 ** L),
                        ess=float(np.exp(-(p * np.log(np.maximum(p, 1e-300))).sum())),
                        top_string=max(strs, key=strs.get),
                        top_frac=float(max(strs.values()) / max(tot, 1)))
    out["ss_window_diversity"] = ssdiv
    print("SS-string diversity: " + json.dumps(
        {L: (v["n_distinct"], round(v["ess"], 1), round(v["top_frac"], 3))
         for L, v in ssdiv.items()}), flush=True)

    out["complete"] = bool(set(out["windows_per_length"]) == set(LENGTHS)
                           and set(ssdiv) == set(LENGTHS)
                           and len(per) >= 20)
    write_json(os.path.join(RESULTS, f"a_charac_cheap_{HASH}.json"), out)
    print(f"cheap done {time.time()-t0:.0f}s")
    return out


# --------------------------------------------------------------------------- heavy
def leaders(W, taus, rng, ngrid):
    """Greedy leader count at each tau, on nested random subsamples.

    Greedy leader clustering: walk the (shuffled) windows; a window becomes a LEADER if its
    Kabsch CA-RMSD to every existing leader exceeds tau.  The leader count is a covering
    number of the set at resolution tau -- the pre-registered 'effective independent count'.

    ADMISSIBLE PREFILTER, so this is exact and not an approximation.  For centred point sets
    a and b, the Kabsch trace obeys  sum(S) <= sqrt(|a|^2 |b|^2), hence
        n * RMSD^2 = |a|^2 + |b|^2 - 2 sum(S) >= (|a| - |b|)^2 = n (Rg_a - Rg_b)^2,
    i.e.  RMSD >= |Rg_a - Rg_b|.  A leader whose radius of gyration differs from the
    candidate's by more than tau therefore CANNOT be within tau, and can be skipped without
    changing the answer.  Verified against the unpruned computation in `_selftest`.
    """
    Wc = W - W.mean(1, keepdims=True)
    rg = np.sqrt((Wc ** 2).sum((1, 2)) / W.shape[1])
    idx = rng.permutation(len(W))
    nmax = min(max(ngrid), len(W))
    out = {}
    for tau in taus:
        lead_i = []
        lead_rg = np.empty(nmax)
        curve = {}
        for c in range(nmax):
            k = idx[c]
            new = True
            if lead_i:
                m = np.abs(lead_rg[:len(lead_i)] - rg[k]) <= tau
                if m.any():
                    sel = np.asarray(lead_i)[m]
                    if I.kabsch_rmsd_batch(Wc[sel], Wc[k]).min() <= tau:
                        new = False
            if new:
                lead_rg[len(lead_i)] = rg[k]
                lead_i.append(k)
            if (c + 1) in ngrid:
                curve[c + 1] = len(lead_i)
        curve[nmax] = len(lead_i)
        out[tau] = curve
    return out


def _selftest():
    """The prefilter must not change the answer.  Brute force vs pruned, small sample."""
    rng = np.random.default_rng(0)
    W = rng.normal(size=(300, 9, 3)) * 3.0
    Wc = W - W.mean(1, keepdims=True)
    for tau in (1.0, 3.0):
        lead = []
        for k in range(len(W)):
            if not lead or I.kabsch_rmsd_batch(Wc[np.array(lead)], Wc[k]).min() > tau:
                lead.append(k)
        got = leaders(W, (tau,), np.random.default_rng(1), (len(W),))
        ref = leaders(W, (tau,), np.random.default_rng(1), (len(W),))
        assert got == ref
        # same walk order as brute force
        brute = len(lead)
        pruned = leaders(W, (tau,), type("R", (), {"permutation": staticmethod(
            lambda n: np.arange(n))})(), (len(W),))[tau][len(W)]
        assert brute == pruned, f"prefilter changed the answer: {brute} vs {pruned}"
    print("prefilter selftest OK", flush=True)


def control_windows(L, n, rng, kind):
    """A PLAUSIBLE zero-information control, in the operator's space (CA coordinates)."""
    if kind == "helix":
        phi = np.full((n, L), math.radians(-57.0))
        psi = np.full((n, L), math.radians(-47.0))
        return I.build_ca(phi, psi)
    if kind == "rama":
        z = np.load(os.path.join(ROOT, "s8", "generate_rama.npz"))
        cnt = np.asarray(z["cnt"], float).sum(0).sum(0)        # marginal over folds and aa
        p = (cnt / cnt.sum()).ravel()
        RB = cnt.shape[0]
        k = rng.choice(len(p), size=n * L, p=p)
        bi, bj = k // RB, k % RB
        step = 2 * math.pi / RB
        phi = (-math.pi + (bi + rng.random(n * L)) * step).reshape(n, L)
        psi = (-math.pi + (bj + rng.random(n * L)) * step).reshape(n, L)
        return I.build_ca(phi, psi)
    raise ValueError(kind)


def heavy():
    fd = os.open(LOCK, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, f"s24 workstream A a_charac heavy pid={os.getpid()}\n".encode())
    os.close(fd)
    try:
        t0 = time.time()
        chains, h = permitted_chains()
        out = dict(CORPUS_HASH=h, taus=list(TAUS), ngrid=list(NGRID), lengths=list(LENGTHS),
                   corpus={}, control_helix={}, control_rama={}, torsion_secondary={})
        for L in LENGTHS:
            w = windows_of(chains, L)
            W = w[0]
            rng = stable_rng("s24A", "leaders", L)
            out["corpus"][L] = dict(n_windows=int(len(W)),
                                    curve=leaders(W, TAUS, rng, NGRID))
            for kind, key in (("helix", "control_helix"), ("rama", "control_rama")):
                C = control_windows(L, max(CTRL_NGRID), stable_rng("s24A", key, L), kind)
                out[key][L] = dict(curve=leaders(C, TAUS, stable_rng("s24A", key, "p", L),
                                                 CTRL_NGRID))
            # SECONDARY (pre-registered): the same clustering in TORSION space, on the
            # same subsample, so a disagreement with the CA-space primary is visible.
            PH, PS = w[1], w[2]
            X = np.concatenate([np.cos(PH), np.sin(PH), np.cos(PS), np.sin(PS)], 1)
            r2 = stable_rng("s24A", "tors", L)
            ii = r2.permutation(len(X))[:max(NGRID)]
            Xs = X[ii]
            lead, tcurve = [], {}
            for c in range(len(Xs)):
                if not lead or np.sqrt(((np.stack(lead) - Xs[c]) ** 2).sum(1)
                                       / (2 * L)).min() > 0.35:
                    lead.append(Xs[c])
                if (c + 1) in NGRID:
                    tcurve[c + 1] = len(lead)
            out["torsion_secondary"][L] = dict(threshold_rad_rms=0.35, curve=tcurve)
            print(f"L={L:2d} n={len(W):7d}  CA leaders@1.0 "
                  f"{out['corpus'][L]['curve'][1.0]}  helix "
                  f"{out['control_helix'][L]['curve'][1.0]}  rama "
                  f"{out['control_rama'][L]['curve'][1.0]}  ({time.time()-t0:.0f}s)",
                  flush=True)
        out["complete"] = bool(set(out["corpus"]) == set(LENGTHS)
                               and all(set(out["corpus"][L]["curve"]) == set(TAUS)
                                       for L in LENGTHS)
                               and set(out["control_rama"]) == set(LENGTHS)
                               and set(out["control_helix"]) == set(LENGTHS)
                               and set(out["torsion_secondary"]) == set(LENGTHS))
        write_json(os.path.join(RESULTS, f"a_charac_heavy_{HASH}.json"), out)
        print(f"heavy done {time.time()-t0:.0f}s")
    finally:
        try:
            os.remove(LOCK)
        except OSError:
            pass


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "cheap"
    if which == "cheap":
        cheap()
    elif which == "heavy":
        heavy()
    else:
        raise SystemExit("cheap|heavy")
