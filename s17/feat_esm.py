"""s17/feat_esm.py -- DO LEARNED SEQUENCE REPRESENTATIONS DISCRIMINATE IN-BAND?

Pre-registration: `s17/PREREG_features.md`, written before any number here.

THE LANE, in one paragraph.  The shipped distance objective has global Spearman 0.568 and IN-BAND
0.131; its argmin inside its own top-25 is -0.014 [-0.118, +0.094], 62W/64L against a random pick
from that same top-25 (LEDGER L10).  Three independent closures -- 58 distogram functionals, 29
native-free in-band signals, 27 target-level calibration features -- all land at zero.  EVERY ONE
of those signals is geometric or energetic.  The one untested feature class is LEARNED SEQUENCE
REPRESENTATION, and the record contains an overturned result claiming ESM buys 0.288 A on
SELECTION (n = 126, p = 0.005) -- measured GLOBALLY, through a retrained distance predictor,
never in-band and never as a per-candidate feature.

THE BAR IS NOT THE DISTANCE OBJECTIVE.  It has no in-band skill, so beating it means nothing.
The bar is: beat RANDOM-IN-BAND, beat a CONSTANT IDEAL ALPHA-HELIX, and beat a CONSTANT IDEAL
BETA-STRAND, on selected CA-RMSD, paired at the target level, with a fold-clustered interval
excluding zero.

THE ARMS.

  F1  ESM-2's CONTACT HEAD against the candidate's realised geometry.  The contact head is the
      only part of ESM-2 supervised on structure (of other proteins; never of the target).  It is
      a DIFFERENT CHANNEL from the shipped distogram, so it is not covered by L5's 58-functional
      closure.
        esmcon_bce    cross-entropy between ESM's contact probability and the candidate's soft
                      realised contact map -- the principled form.
        esmcon_agree  contact-probability-weighted mean realised contact.
        esmcon_corr   -correlation between the two maps.
        esmcon_topd   mean realised CA distance over the n pairs ESM ranks most contact-like.

  F2  ESM contact CONFIDENCE as an EXOGENOUS per-pair WEIGHT on the shipped Bayes risk.  L5
      closed the distogram's own weightings (sd, entropy, sequence separation).  A weight derived
      from a different model is the one thing that sweep did not contain.
        risk_wcon / risk_wconf / risk_wshipcon, against risk_wone (uniform) and dist_shipped.

  F4  DOES THE EMBEDDING PREDICT THIS RESIDUE'S LOCAL ENVIRONMENT?  Leave-fold-out ridge from the
      per-residue ESM-PCA32 block (with +-1 neighbours and position) to the residue's native local
      CA descriptor [d(r,r+2), d(r,r+3), d(r,r+4)].  Candidate score = standardised disagreement
      with its realised descriptor.
        env_esm   vs   env_onehot  (THE CONTROL THAT DECIDES IT -- the exact ESM/one-hot contrast
                                    the 0.288 A claim rests on)
                  vs   env_const   (training mean, no sequence at all: zero-information)

  F5  PREDICTED vs REALISED SECONDARY STRUCTURE, per residue.  Same three arms, multinomial
      logistic, 3-class CA-only SS assigned by ONE rule applied identically to native and
      candidate.
        ss_esm  vs  ss_onehot  vs  ss_marg

  F3  Pseudo-likelihood of the candidate window's OWN source sequence under the target's ESM
      context (`feat_pll.py`).  Included automatically when its cache exists.

LEAKAGE.  F4/F5 labels are native structures of TRAINING-FOLD targets only; the held-out target's
native never enters its own score.  Asserted by construction: the fold loop fits on `fold != f`
and predicts on `fold == f`.  `rr` appears only in evaluation.

MIN-OF-N.  Every feature's sign is fixed a priori, so this battery is not a hidden best-of-2V
screen.  The report nevertheless prints the L8 zero-information min-of-N null beside any
per-target best-of-V statistic, and quotes no such statistic without it.
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

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s17 import sel_lib as L               # noqa: E402
from s17 import feat_lib as F              # noqa: E402

K = 500
OUTJSON = os.path.join(F.RESULTS, "feat_inband.json")
#: which per-candidate window-embedding cache F6 consumes.  Declared here rather than discovered,
#: so the artefact records which instrument produced the numbers (a 650M model and an 8M model are
#: different instruments and their results must never be pooled).
FW_B = int(os.environ.get("S17_FW_B", "25"))
FW_MODEL = os.environ.get("S17_FW_MODEL", "esm2_t33_650M_UR50D")
BAND_ONLY = ("emb_cos", "emb_l2", "emb_pool")


# --------------------------------------------------------------------------- per-residue blocks
def _block(X, n):
    """(n, 3d + 3) context block: residue r with its +-1 neighbours, plus position features."""
    d = X.shape[1]
    pad = np.zeros((1, d), X.dtype)
    Xm = np.vstack([pad, X[:-1]])
    Xp = np.vstack([X[1:], pad])
    pos = np.column_stack([np.arange(n) / max(n - 1, 1),
                           (n - 1 - np.arange(n)) / max(n - 1, 1),
                           np.full(n, n / 16.0)])
    return np.column_stack([Xm, X, Xp, pos]).astype(np.float64)


def _rows_for(t, arm):
    """(m, p) design block and (m, 3) native descriptor for the first n-4 residues."""
    seq, n = t["seq"], int(t["n"])
    X = F.esm_pca32(seq) if arm == "esm" else F.onehot(seq)
    B = _block(np.asarray(X, np.float64), n)
    u = I.load_univ(t["pdb"])
    nat = np.asarray(u["nat_ca"], float)
    E = F.env_desc(nat)                      # (n-4, 3)
    S = F.ss_ca(nat)                         # (n,)
    m = len(E)
    return B[:m], E, S[:m]


def _ridge(A, Y, lam=10.0):
    A1 = np.column_stack([A, np.ones(len(A))])
    G = A1.T @ A1 + lam * np.eye(A1.shape[1])
    G[-1, -1] -= lam                          # do not penalise the intercept
    return np.linalg.solve(G, A1.T @ Y)


def _apply(A, w):
    return np.column_stack([A, np.ones(len(A))]) @ w


def fit_heads(tg, verbose=True):
    """Leave-fold-out per-residue heads.  Labels are natives of TRAINING folds only."""
    from sklearn.linear_model import LogisticRegression
    cache = {}
    for arm in ("esm", "onehot"):
        for t in tg:
            cache[(arm, t["pdb"])] = _rows_for(t, arm)
    folds = np.array([int(t["fold"]) for t in tg], int)
    heads = {}
    for f in sorted(set(folds.tolist())):
        tr = [t for t, g in zip(tg, folds) if g != f]
        for arm in ("esm", "onehot"):
            A = np.vstack([cache[(arm, t["pdb"])][0] for t in tr])
            Y = np.vstack([cache[(arm, t["pdb"])][1] for t in tr])
            S = np.concatenate([cache[(arm, t["pdb"])][2] for t in tr])
            mu, sd = Y.mean(0), Y.std(0) + 1e-9
            heads[(f, arm, "env")] = (_ridge(A, (Y - mu) / sd), mu, sd)
            lr = LogisticRegression(C=0.5, max_iter=400)
            classes = np.unique(S)
            if len(classes) > 1:
                lr.fit(A, S)
                heads[(f, arm, "ss")] = lr
            else:
                heads[(f, arm, "ss")] = None
        #: ZERO-INFORMATION arms: the training mean descriptor and the training SS marginal,
        #: with no sequence input whatsoever.
        Y = np.vstack([cache[("esm", t["pdb"])][1] for t in tr])
        S = np.concatenate([cache[("esm", t["pdb"])][2] for t in tr])
        heads[(f, "const", "env")] = (Y.mean(0), Y.mean(0), Y.std(0) + 1e-9)
        marg = np.array([(S == c).mean() for c in range(3)])
        heads[(f, "const", "ss")] = np.log(np.maximum(marg, 1e-6))
        if verbose:
            print(f"  heads fold {f}: {len(tr)} training targets, "
                  f"{sum(len(cache[('esm', t['pdb'])][1]) for t in tr)} residues", flush=True)
    return heads


# --------------------------------------------------------------------------- per-candidate features
def features(t, p, heads, sc):
    """Every representation-derived per-candidate feature.  LOWER = PREDICTED BETTER, always."""
    seq, n, f = t["seq"], int(t["n"]), int(t["fold"])
    D = np.asarray(p["D"], float)
    W = np.asarray(p["W"], float)
    i, j = np.asarray(p["i"], int), np.asarray(p["j"], int)
    out = {}

    # ---------------------------------------------------------------- F1  ESM contact head
    con = F.esm_contacts()[seq]
    pc = np.clip(np.asarray(con, float)[i, j], 1e-6, 1 - 1e-6)
    s = np.clip(F.soft_contact(D), 1e-6, 1 - 1e-6)
    out["esmcon_bce"] = (-(pc[None, :] * np.log(s) + (1 - pc[None, :]) * np.log(1 - s))).mean(1)
    out["esmcon_agree"] = -(s * pc[None, :]).sum(1) / max(pc.sum(), 1e-9)
    sz = (s - s.mean(1, keepdims=True)) / (s.std(1, keepdims=True) + 1e-9)
    pz = (pc - pc.mean()) / (pc.std() + 1e-9)
    out["esmcon_corr"] = -(sz * pz[None, :]).mean(1)
    top = np.argsort(-pc)[:max(n, 1)]
    out["esmcon_topd"] = D[:, top].mean(1)

    #: THE CONTROL THAT DECIDES WHETHER F1 IS ESM AT ALL.  `esmcon_agree` is a weighted count of
    #: realised contacts, and a weighted contact count is dominated by HOW COMPACT the candidate
    #: is.  Its ESM-free twin replaces p_esm by a uniform weight; `rg_raw` is the crudest
    #: compactness statistic there is.  If the ESM-weighted form does not beat these, F1 is
    #: measuring compactness -- which SELECT already surveyed -- and not the language model.
    #: This is BRIEF s5's failure mode (a quantity measured correctly, read as a different
    #: quantity) and it is the one place in this lane where it could have landed.
    out["cf_uniform"] = -s.mean(1)
    #: matched to `esmcon_topd`: the n pairs the DISTOGRAM thinks are most contact-like, i.e. the
    #: smallest predicted distances, rather than the n pairs ESM's contact head ranks highest.
    out["cf_topd_uniform"] = D[:, np.argsort(p["expected"])[:max(n, 1)]].mean(1)
    out["rg_raw"] = np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(-1).mean(-1))

    # ---------------------------------------------------------------- F2  exogenous pair weights
    R = p["risk_raw"][np.arange(p["risk_raw"].shape[0])[None, :], L.gbin(D, p["grid"])]
    for nm, w in (("risk_wone", np.ones_like(pc)),
                  ("risk_wcon", pc),
                  ("risk_wconf", np.maximum(pc, 1 - pc)),
                  ("risk_wshipcon", p["w_ship"] * pc)):
        wn = w / max(np.mean(w), 1e-12)
        out[nm] = (R * wn[None, :]).mean(1)
    out["dist_shipped"] = sc                     # REFERENCE, not a bar

    # ---------------------------------------------------------------- F4  local environment
    Ec = F.env_desc(W)                           # (B, n-4, 3)
    m = Ec.shape[1]
    for arm in ("esm", "onehot", "const"):
        w, mu, sd = heads[(f, arm, "env")]
        if arm == "const":
            pred = np.tile(((mu - mu) / sd), (m, 1))
        else:
            X = F.esm_pca32(seq) if arm == "esm" else F.onehot(seq)
            pred = _apply(_block(np.asarray(X, np.float64), n)[:m], w)
        Z = (Ec - mu[None, None, :]) / sd[None, None, :]
        out[f"env_{arm}"] = np.abs(Z - pred[None, :, :]).mean((1, 2))

    # ---------------------------------------------------------------- F5  secondary structure
    Sc = F.ss_ca(W)[:, :m]                       # (B, n-4)
    for arm in ("esm", "onehot", "const"):
        h = heads[(f, arm, "ss")]
        if arm == "const":
            lp = np.tile(h[None, :], (m, 1))
        elif h is None:
            lp = np.zeros((m, 3))
        else:
            X = F.esm_pca32(seq) if arm == "esm" else F.onehot(seq)
            P = h.predict_proba(_block(np.asarray(X, np.float64), n)[:m])
            lp = np.full((m, 3), np.log(1e-6))
            for q, c in enumerate(h.classes_):
                lp[:, int(c)] = np.log(np.maximum(P[:, q], 1e-6))
        out[f"ss_{arm}"] = -lp[np.arange(m)[None, :], Sc].mean(1)

    # ---------------------------------------------------------------- F3  window-sequence PLL
    #: the SEQUENCE-CHANNEL CONTROLS: BLOSUM62 similarity (a position-independent substitution
    #: matrix -- the classical answer to "does this fragment's sequence match") and the
    #: retrieval rank it induces.  Both already measured null in-band by SELECT; carried here so
    #: the PLL arm has a paired contrast rather than only an absolute number.
    u_ = I.load_univ(t["pdb"])
    out["blosum_sim"] = -np.asarray(u_["sim"], float)[p["idx"]]
    from s17 import feat_pll as PL
    pll = PL.load().get(seq)
    if pll is not None:
        Scodes = np.asarray(u_["S"], np.int64)[p["idx"]]                      # (B, n)
        for nm, lp in (("pll_nat", pll[0]), ("pll_mask", pll[1])):
            out[nm] = -lp[np.arange(n)[None, :], Scodes].mean(1)
        #: ZERO-INFORMATION control: the same score under a position-INDEPENDENT background
        #: frequency, which knows the amino-acid composition and nothing else.
        bg = np.bincount(Scodes.ravel(), minlength=20).astype(float)
        bg = np.log(np.maximum(bg / bg.sum(), 1e-6))
        out["pll_bg"] = -bg[Scodes].mean(1)

    # ---------------------------------------------------------------- F6  window embeddings
    #: BAND-RESTRICTED.  Per-candidate ESM embeddings exist only for the band `feat_window.py`
    #: featurised (a cost gate declared in the pre-registration), so these features are valid on
    #: exactly those candidates and are reported for no other band.  `valid` carries that.
    from s17 import feat_window as FW
    valid = None
    wc = FW.load(B=FW_B, model=FW_MODEL)
    if wc is not None and t["pdb"] in wc["band"]:
        rk, wseqs = wc["band"][t["pdb"]]
        valid = np.zeros(len(sc), bool)
        valid[rk] = True
        Et = np.asarray(F.esm_reps()[seq], np.float64)                   # (n, 1280)
        Zt = np.asarray(F.esm_pca32(seq), np.float64)                    # (n, 32)
        mu, Wp, scale = F._PCA
        big = np.full(len(sc), np.nan)
        cos = np.full(len(sc), np.nan)
        pool = np.full(len(sc), np.nan)
        tn = Et / np.maximum(np.linalg.norm(Et, axis=1, keepdims=True), 1e-9)
        tp = Et.mean(0)
        tp = tp / max(np.linalg.norm(tp), 1e-9)
        for q, (ci, ws) in enumerate(zip(rk, wseqs)):
            raw = wc["emb"].get(ws)
            if raw is None:
                continue
            Ew = np.asarray(raw, np.float64)
            if Ew.ndim != 2 or Ew.shape[0] != n:
                continue
            wn = Ew / np.maximum(np.linalg.norm(Ew, axis=1, keepdims=True), 1e-9)
            cos[ci] = -float((tn * wn).sum(1).mean())
            #: the whitened-PCA distance is only defined when the window embedding lives in the
            #: SAME space the projection was fitted in.  A smaller scout model (8M, 320-d) does
            #: not, so `emb_l2` is simply absent there rather than silently computed in a space
            #: it does not belong to -- the cosine arms are dimension-agnostic and still run.
            if Ew.shape[1] == len(mu):
                Zw = ((Ew - mu) @ Wp) / np.maximum(scale, 1e-6)
                big[ci] = float(np.linalg.norm(Zt - Zw, axis=1).mean())
            wp = Ew.mean(0)
            pool[ci] = -float(tp @ (wp / max(np.linalg.norm(wp), 1e-9)))
        out["emb_cos"] = cos
        out["emb_pool"] = pool
        if not np.isnan(big).all():
            out["emb_l2"] = big
        valid &= ~np.isnan(cos)
    return out, valid


def one_target(t, heads):
    p = L.pack(t["pdb"], K, want=("D", "W"))
    sc = L.shipped(p)
    rr = p["rr"]
    W = p["W"]
    rng = SD.stable_rng(t["pdb"], "s17feat")
    sg, valid = features(t, p, heads, sc)
    sg.update(F.null_signals(p, W, rng))
    masks = F.bands_of(sc, rr)
    row = {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]),
           "oracle": float(rr.min()), "global": {}, "band": {}}
    for nm, v in sg.items():
        if nm not in BAND_ONLY:
            row["global"][nm] = F.spearman(v, rr)
    for b, mk in masks.items():
        cell = {"size": int(mk.sum()), "best": float(rr[mk].min()),
                "mean": float(rr[mk].mean()), "sig": {}}
        if mk.sum() >= 8:
            for nm, v in sg.items():
                #: a band-restricted feature is scored only where every band member has it --
                #: scoring it on a partial band would compare it against a DIFFERENT candidate
                #: set from every other feature, which is the one thing the instrument forbids.
                if nm in BAND_ONLY:
                    if valid is None or not valid[mk].all():
                        continue
                cell["sig"][nm] = F.inband_stats(v[mk], rr[mk])
            #: a LADDER of V, not a single value: the number of features present varies with which
            #: caches exist, and the report must be able to match V to the feature count it
            #: actually compares without the artefact having to be regenerated.
            cell["minN"] = {str(V): F.minN_null(rr[mk], V, t["pdb"], f"{b}_{V}")
                            for V in (2, 3, 5, 8, 10, 15, 19, 20, 22, 25, 30)}
        row["band"][b] = cell
    return row


def run(targets=None, verbose=True):
    tg = targets if targets is not None else I.targets()
    print(f"fitting leave-fold-out heads on {len(tg)} targets ...", flush=True)
    heads = fit_heads(tg, verbose=verbose)
    rows, t0 = [], time.time()
    for q, t in enumerate(tg):
        rows.append(one_target(t, heads))
        if verbose and (q + 1) % 10 == 0:
            print(f"  {q+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"rows": rows}, open(OUTJSON, "w"))
    json.dump({"rows": rows}, open(OUTJSON, "w"))
    print(f"wrote {OUTJSON}  ({time.time()-t0:.0f}s)", flush=True)
    return rows


if __name__ == "__main__":
    run()
