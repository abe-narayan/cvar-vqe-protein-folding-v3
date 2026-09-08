"""SPRINT 20, AGENT A -- CAN ANY CANDIDATE GENERATOR DECORRELATE FROM THE RETRIEVAL POOL?

Pre-registration: `s20/PREREG_A.md`, written before this module ran and not edited since.

THE OPERATOR (identical for every arm; only the candidate index set moves):

    rank a window subset by the universe's own stable BLOSUM `order` (filtered, so ties break
    identically)  ->  first 500  ->  shipped Bayes-risk distogram score  ->  lowest 75
    ->  X_S = I.coordinate_average(top-75)                  the EMITTED structure (raw)

BASIS -- CORRECTED 2026-09-07 by the coordinator, against this module's own pre-registration.
`I.coordinate_average` returns a POINT CLOUD, not a buildable chain: its mean virtual Ca-Ca bond is
2.961 A against the physical 3.804 A (22.2% contraction, global min bond 0.649 A).  My
pre-registration (`s20/PREREG_A.md` section 0) called it "the best built object" and made it the
primary emitted object.  That was wrong and the pre-registration is left unedited with the
deviation recorded in `s20/agentA_FINDINGS.md`.

  * every point-cloud number here is labelled `rmsd|` and is comparable ONLY to other point clouds
    (reference: coordinate average of the shipped top-75, 3.048 A)
  * the BUILT-CHAIN number is `a_proj.py`'s `rmsdP|`, the lam=0 projection, comparable to the
    incumbent 3.204 A
  * `bond|` / `bondmin|` / `bondmem|` carry the virtual Ca-Ca bond so no arm can buy RMSD by
    contracting the chain

THE ERROR FIELD.  X_S is a real structure, so

    e_S = d(X_S) - d_true          (min_sep = 2)

IS its coherent error: the Sprint-19 split r = r_coh + r_inc is degenerate here with r_inc = 0.
That is an IDENTITY, not a finding.  rho(S,T) is the mean-centred Pearson correlation of e_S and
e_T across a target's pairs, averaged over targets -- the same semantics as `s19/a_source._corr`.

THE MATCHED CONTROL, in the operator's own space.  `pep`/`prot` is a partition of the universe by
CORPUS.  `halfA`/`halfB` is a partition of the same universe AT RANDOM, consumed by the same
operator at the same counts.  The corpus can only be the carrier of the shared bias if the corpus
partition decorrelates MORE than the random partition.

ORACLE: `dtrue`, `rr`, `nat_ca` and everything derived from them.  Evaluation only; never enters a
candidate set or a score.

Run:  python -m s20.a_src
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

from s12 import instrument as I              # noqa: E402
from s15 import seed as SD                   # noqa: E402

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(HERE, "cache")
for _d in (RESULTS, CACHE):
    os.makedirs(_d, exist_ok=True)

OUT = os.path.join(RESULTS, "a_src.json")

K = 500
M = 75

RETRIEVAL = ["pool", "pep", "prot", "halfA", "halfB", "rand500"]
GENERATED = ["tors", "unsel"]
SOURCES = RETRIEVAL + GENERATED
FUSE = [("pool", "pep"), ("pool", "prot"), ("pool", "tors"), ("pool", "rand500"),
        ("pool", "unsel"), ("halfA", "halfB")]
LONG_SEP = 5          # the sep >= 5 restriction for the long-range variant of rho

# ---- classical continuous-torsion generator: 3-state first-order Markov chain, per-state
# Gaussian (phi, psi).  Parameters are standard backbone geometry, NOT fitted to this project's
# corpus, and the only input is the chain length.  Sequence-blind, zero target information.
_MU = np.deg2rad(np.array([[-63.0, -42.0],       # H  alpha
                           [-125.0, 135.0],      # E  beta
                           [-80.0, 130.0]]))     # C  coil / PPII
_SG = np.deg2rad(np.array([[8.0, 8.0], [25.0, 25.0], [35.0, 40.0]]))
_TRANS = np.array([[0.88, 0.02, 0.10],
                   [0.02, 0.85, 0.13],
                   [0.10, 0.13, 0.77]])
_INIT = np.array([0.35, 0.20, 0.45])
_LMU = np.deg2rad(np.array([57.0, 40.0]))        # left-handed alpha, 8% of coil residues
_LSG = np.deg2rad(np.array([12.0, 12.0]))


def torsion_sample(n, b, rng):
    """(b, n) phi and psi from the sequence-blind Markov torsion model."""
    st = np.empty((b, n), np.int64)
    st[:, 0] = rng.choice(3, size=b, p=_INIT)
    for k in range(1, n):
        u = rng.random(b)
        cum = np.cumsum(_TRANS[st[:, k - 1]], axis=1)
        st[:, k] = (u[:, None] > cum).sum(1)
    mu = _MU[st]
    sg = _SG[st]
    x = mu + sg * rng.standard_normal((b, n, 2))
    left = (st == 2) & (rng.random((b, n)) < 0.08)
    xl = _LMU + _LSG * rng.standard_normal((b, n, 2))
    x = np.where(left[..., None], xl, x)
    return x[..., 0], x[..., 1]


def _score(dg, W, i, j):
    D = I.pair_dists(W, i, j)
    return I.shipped_score(dg, D.astype(np.float32).astype(float)), D


def _std(E):
    """Row-standardise an (m, npairs) error matrix so E @ E.T is the correlation matrix."""
    E = np.asarray(E, float)
    E = E - E.mean(1, keepdims=True)
    nrm = np.linalg.norm(E, axis=1, keepdims=True)
    return E / np.maximum(nrm, 1e-12)


def run(tg=None, limit=None):
    tg = tg if tg is not None else I.targets()
    if limit:
        tg = tg[:limit]
    rows = []
    gate_ok = 0
    t0 = time.time()
    for c, t in enumerate(tg):
        pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
        u = I.load_univ(pdb)
        nat = u["nat_ca"]
        dg = I.distogram(pdb, seq, fold)
        i, j = I.pair_index(n)
        dtrue = np.linalg.norm(nat[i] - nat[j], axis=1)
        order = np.asarray(u["order"], int)
        org = np.asarray(u["org"], bool)
        nw = len(org)
        W_all = np.asarray(u["W"], float)
        rr = np.asarray(u["rr"], float)

        rng = SD.stable_rng(pdb, "s20A_half")
        half = rng.random(nw) < 0.5

        idx = {"pool": order[:K],
               "pep": order[org[order]][:K],
               "prot": order[~org[order]][:K],
               "halfA": order[half[order]][:K],
               "halfB": order[~half[order]][:K],
               "rand500": np.sort(SD.stable_rng(pdb, "s20A_rand")
                                  .choice(nw, size=min(K, nw), replace=False))}

        e = {"pdb": pdb, "n": n, "fold": fold, "npairs": int(len(i)),
             "nw": int(nw), "pep_frac_univ": float(org.mean()),
             "pep_frac_pool": float(org[order[:K]].mean())}

        X, TOP, EMEM, POOLRR = {}, {}, {}, {}
        for s in RETRIEVAL:
            ix = idx[s]
            W = W_all[ix]
            sc, D = _score(dg, W, i, j)
            sub = np.argsort(sc, kind="stable")[:M]
            top = ix[sub]
            if s == "pool":
                gate_ok += int(set(top.tolist()) ==
                               set(order[:K][np.asarray(I.shipped_record(pdb)["sub"], int)].tolist()))
            Wt = W_all[top]
            Xs, _b = I.coordinate_average(Wt)
            X[s] = Xs
            TOP[s] = top
            EMEM[s] = D[sub] - dtrue[None]          # ORACLE member coherent error fields
            POOLRR[s] = rr[ix]
            e[f"npool|{s}"] = int(len(ix))
            e[f"gen_ceiling|{s}"] = float(rr[ix].min())
            e[f"sel_ceiling|{s}"] = float(rr[top].min())
            e[f"member_rmsd|{s}"] = float(rr[top].mean())
            e[f"cov2|{s}"] = float((rr[top] < 2.0).mean())
            e[f"cov2pool|{s}"] = float((rr[ix] < 2.0).mean())

        # ---- generated source: the classical torsion model
        rngt = SD.stable_rng(pdb, "s20A_tors")
        PHI, PSI = torsion_sample(n, K, rngt)
        Wg = np.asarray(I.build_ca(PHI, PSI), float)
        sc, D = _score(dg, Wg, i, j)
        sub_t = np.argsort(sc, kind="stable")[:M]
        sub = sub_t
        Wt = Wg[sub]
        rrg = I.kabsch_rmsd_batch(Wg, nat)
        X["tors"], _b = I.coordinate_average(Wt)
        TOP["tors"] = None
        EMEM["tors"] = D[sub] - dtrue[None]
        e["npool|tors"] = int(K)
        e["gen_ceiling|tors"] = float(rrg.min())
        e["sel_ceiling|tors"] = float(rrg[sub].min())
        e["member_rmsd|tors"] = float(rrg[sub].mean())
        e["cov2|tors"] = float((rrg[sub] < 2.0).mean())
        e["cov2pool|tors"] = float((rrg < 2.0).mean())

        # ---- SCORE-FREE ZERO-INFORMATION RETRIEVAL REFERENCE: 75 uniformly random windows
        # from the WHOLE universe -- no BLOSUM, no distogram score.  This is the reference for
        # "how much do two n-residue peptides that are not this native agree on being wrong".
        iu75 = SD.stable_rng(pdb, "s20A_unsel").choice(nw, size=M, replace=False)
        Wu75 = W_all[iu75]
        X["unsel"], _b = I.coordinate_average(Wu75)
        TOP["unsel"] = iu75
        EMEM["unsel"] = I.pair_dists(Wu75, i, j) - dtrue[None]
        e["npool|unsel"] = int(nw)
        e["gen_ceiling|unsel"] = float(rr.min())
        e["sel_ceiling|unsel"] = float(rr[iu75].min())
        e["member_rmsd|unsel"] = float(rr[iu75].mean())
        e["cov2|unsel"] = float((rr[iu75] < 2.0).mean())
        e["cov2pool|unsel"] = float((rr < 2.0).mean())

        WT = {s: (W_all[TOP[s]] if TOP[s] is not None else Wt) for s in SOURCES}

        # ---- zero-information reference structure
        X["helix"] = I.build_ca(np.full(n, np.deg2rad(-63.0)), np.full(n, np.deg2rad(-42.0)))

        # ---- emitted structures, their errors, diversity, and the VIRTUAL BOND axis
        E = {}
        for s in list(SOURCES) + ["helix"]:
            e[f"rmsd|{s}"] = float(I.ca_rmsd(X[s], nat))
            E[s] = np.linalg.norm(X[s][i] - X[s][j], axis=1) - dtrue
            e[f"err_rms|{s}"] = float(np.sqrt((E[s] ** 2).mean()))
            bl = np.linalg.norm(X[s][1:] - X[s][:-1], axis=1)
            e[f"bond|{s}"] = float(bl.mean())
            e[f"bondmin|{s}"] = float(bl.min())
        e["bond|native"] = float(np.linalg.norm(nat[1:] - nat[:-1], axis=1).mean())
        for s in SOURCES:
            blm = np.linalg.norm(WT[s][:, 1:] - WT[s][:, :-1], axis=2)
            e[f"bondmem|{s}"] = float(blm.mean())
            e[f"bondmemmin|{s}"] = float(blm.min())
            P = I.pairwise_rmsd(WT[s])
            iu = np.triu_indices(len(P), 1)
            e[f"div_geo|{s}"] = float(P[iu].mean())
            Z = _std(EMEM[s])
            Cm = Z @ Z.T
            e[f"div_err|{s}"] = float(1.0 - Cm[iu].mean())
            e[f"mem_err_rms|{s}"] = float(np.sqrt((EMEM[s] ** 2).mean()))

        # ---- rho(e_coherent, e_pool), the full cross matrix, and the cross-source RMSD
        names = list(SOURCES) + ["helix"]
        Es = _std(np.stack([E[s] for s in names]))
        Cx = Es @ Es.T
        lg = (j - i) >= LONG_SEP
        El = _std(np.stack([E[s][lg] for s in names])) if lg.sum() >= 4 else None
        Cl = El @ El.T if El is not None else None
        for a in range(len(names)):
            for b in range(a + 1, len(names)):
                e[f"rho|{names[a]}|{names[b]}"] = float(Cx[a, b])
                e[f"rhoL|{names[a]}|{names[b]}"] = float(Cl[a, b]) if Cl is not None else np.nan
                e[f"xrmsd|{names[a]}|{names[b]}"] = float(I.ca_rmsd(X[names[a]], X[names[b]]))
        e["n_long"] = int(lg.sum())

        # ---- cross-source MEMBER error correlation (error diversity ACROSS populations)
        for a, b in [("pep", "prot"), ("halfA", "halfB"), ("pool", "pep"), ("pool", "prot"),
                     ("pool", "tors"), ("pool", "rand500"), ("pool", "unsel"),
                     ("unsel", "unsel2")]:
            if b == "unsel2":
                continue
            Za, Zb = _std(EMEM[a]), _std(EMEM[b])
            e[f"memrho|{a}|{b}"] = float((Za @ Zb.T).mean())
        for s in SOURCES:
            Z = _std(EMEM[s])
            iu = np.triu_indices(len(Z), 1)
            e[f"memrho|{s}|{s}"] = float((Z @ Z.T)[iu].mean())

        # ---- fusion arms (POINT CLOUD basis)
        for a, b in FUSE:
            xb = I.superpose_batch(X[b][None], X[a])[0]
            xf = 0.5 * (X[a] + xb)
            e[f"fuse|{a}|{b}"] = float(I.ca_rmsd(xf, nat))
            e[f"fusebond|{a}|{b}"] = float(np.linalg.norm(xf[1:] - xf[:-1], axis=1).mean())
            Wu = np.concatenate([WT[a], WT[b]], 0)
            xm, _ = I.coordinate_average(Wu)
            e[f"merge|{a}|{b}"] = float(I.ca_rmsd(xm, nat))
            e[f"mergebond|{a}|{b}"] = float(np.linalg.norm(xm[1:] - xm[:-1], axis=1).mean())

        # ---- persist the candidate sets so the BUILT-CHAIN pass (a_proj) needs no rescoring
        np.savez(os.path.join(CACHE, f"sets_{pdb}.npz"),
                 **{f"top_{s}": np.asarray(TOP[s] if TOP[s] is not None else [], int)
                    for s in SOURCES},
                 tors_phi=PHI[sub_t], tors_psi=PSI[sub_t])

        rows.append(e)
        if (c + 1) % 10 == 0 or c + 1 == len(tg):
            print(f"  {c+1}/{len(tg)}  ({time.time()-t0:.0f}s)  gate {gate_ok}/{c+1}", flush=True)
            json.dump({"rows": rows, "gate_ok": gate_ok, "complete": len(rows) == len(tg)},
                      open(OUT, "w"), default=float)
    json.dump({"rows": rows, "gate_ok": gate_ok, "complete": len(rows) == len(tg)},
              open(OUT, "w"), default=float)
    if len(rows) == len(tg) and gate_ok == len(tg):
        open(os.path.join(RESULTS, "a_src.COMPLETE"), "w").write("ok\n")
    return rows


if __name__ == "__main__":
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else None
    run(limit=lim)
