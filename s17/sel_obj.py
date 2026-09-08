"""s17/sel_obj.py -- TAKING THE DISTANCE OBJECTIVE APART.  Sprint 17 section 12, PREREG E2.

The distance objective is the only selection signal in this programme that has ever beaten
random with an interval excluding zero (certified argmin -0.926 A [-1.440, -0.431] on 19
enumerated spaces, where neither Legacy nor AMBER is distinguishable from random).  It is
therefore worth knowing what it actually is, and whether the shipped functional is the right
one.

WHAT THE SHIPPED OBJECTIVE IS, read out of `core/predict.py` rather than assumed:

    score(x) = mean_p  w_p * Risk_p(d_p(x)),      Risk_p(d) = sum_b prob_pb * |d - centre_b|
    w_p      = 1 / (sd_p + 0.5),  normalised to mean 1          (SHELLS all weight 1.0,
                                                                 gamma = 1.0 -- checked)

So it is an **L1 Bayes risk**, **uncertainty-weighted**, **unweighted in sequence separation**,
and aggregated by an **unweighted arithmetic mean over pairs**.  Four choices, none of which
has been priced against its alternatives on a fixed candidate set.

THE VARIANT SPACE.  A variant is (loss form) x (pair weight) x (aggregation), plus three
families that do not factor that way:

  loss     risk | abs | z | huber(z) | tukey(z) | log | nll | sq
  weight   one | ship | ship^2 | sqrt(ship) | entropy | sep | 1/sep | long-range | short-range
  agg      mean | median | trim25 | cvar25 | zset | rankset

  centre-blend   score_beta with the predicted centre replaced by the candidate SET's own
                 centre in a one-parameter blend; beta = 0 is the distance objective and
                 beta = 1 is TYPICALITY, the only signal in the record with positive in-band
                 skill.  The interior has never been measured.
  contact        the distance objective collapsed to a contact/no-contact likelihood
  spectral       the low-rank / Torgerson form: match the eigenspectrum of the double-centred
                 squared-distance matrix rather than the distances themselves

TRIANGLE-INEQUALITY AND CHAIN-CONNECTIVITY CONSISTENCY TERMS ARE DELIBERATELY ABSENT, and
that is a result rather than an omission: every candidate here is a real protein backbone
window, so it satisfies the triangle inequality and ideal CA-CA connectivity exactly by
construction.  Such a term is identically zero on this candidate set and cannot rank it.
It would only matter for a generator that emits distance matrices directly.

METHODOLOGY THAT IS NOT NEGOTIABLE HERE.

* Every variant sees the IDENTICAL candidate set.
* `sel_lib.sel_of` averages over the tied argmin set, so a degenerate score scores the set
  mean rather than reading the BLOSUM order.
* The headline is the LEAVE-FOLD-OUT variant choice.  The in-sample best-of-N is printed
  beside it purely as the overfitting upper bound.  `score-axis-does-not-transfer` records
  eleven variants that were promoted on a best-of-N screen and died; that screen had no
  fold-honest selection step and this one does.
* Every variant reports rho as well as RMSD, because rho is what the sprint's design
  equation consumes: `rho_within` (mean per-target Spearman of score against true RMSD --
  the copula's parameter), `rho_band` (the same restricted to the near-native band, the
  quantity comparable with the recorded 0.600), and `rho_pool` (pooled across targets, which
  measures whether the score is comparable BETWEEN targets at all).

PRE-REGISTERED FALSIFIER (PREREG_select.md E2).  If the leave-fold-out variant's advantage
over the shipped objective has a CI containing zero, the space of functionals of the existing
distogram is closed for SELECTION as well as for generation, and that null is the headline.
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

from s17 import sel_lib as L               # noqa: E402

EPS = 1e-9


# ------------------------------------------------------------------ per-target context
def context(p):
    """Everything every variant needs, computed once per target."""
    D = np.asarray(p["D"], float)                       # (K, P)
    K, P = D.shape
    exp = p["expected"]; sd = p["sd"]; sep = p["sep"]
    g = L.gbin(D, p["grid"])
    rows = np.arange(P)[None, :]

    c = {}
    c["D"] = D
    c["risk"] = p["risk_raw"][rows, g]                  # (K, P) unweighted Bayes risk
    r = D - exp[None, :]
    c["abs"] = np.abs(r)
    c["sq"] = r * r
    c["z"] = np.abs(r) / np.maximum(sd, 1e-3)[None, :]
    c["log"] = np.abs(np.log(np.maximum(D, 1e-3)) - np.log(np.maximum(exp, 1e-3))[None, :])
    # Huber and Tukey on the z-residual: bounded influence of a badly-predicted pair
    z = c["z"]
    c["hub"] = np.where(z <= 1.0, 0.5 * z * z, z - 0.5)
    zc = np.minimum(z, 3.0) / 3.0
    c["tuk"] = 1.0 - (1.0 - zc * zc) ** 3
    # NLL of the observed bin under the predicted distribution
    edges = 0.5 * (p["centres"][1:] + p["centres"][:-1])
    b = np.digitize(D, edges)
    c["nll"] = -np.log(np.maximum(p["prob"][rows, b], 1e-12))

    c["setmean"] = D.mean(0)
    c["setsd"] = np.maximum(D.std(0), 1e-3)
    c["exp"] = exp; c["sd"] = sd; c["sep"] = sep; c["K"] = K; c["P"] = P
    c["prob"] = p["centres"], p["prob"]

    # pair weights, each normalised to mean 1 so aggregation scales are comparable
    ent = -(p["prob"] * np.log(np.maximum(p["prob"], 1e-12))).sum(1)
    ship = p["w_ship"]
    def nz(w):
        w = np.asarray(w, float)
        return w / max(w.mean(), 1e-12)
    c["w"] = {
        "one": np.ones(P),
        "ship": ship,
        "ship2": nz(ship ** 2),
        "shiphalf": nz(np.sqrt(ship)),
        "ent": nz(1.0 / (ent + 0.5)),
        "maxp": nz(p["prob"].max(1)),
        "sep": nz(sep),
        "isep": nz(1.0 / sep),
        "lr": nz((sep >= 6).astype(float) + 1e-6),
        "sr": nz((sep <= 5).astype(float) + 1e-6),
        "shipsep": nz(ship * sep),
        "shipisep": nz(ship / sep),
    }
    return c


# ------------------------------------------------------------------ aggregations
def _agg(Lm, w, how):
    """Aggregate a (K, P) per-pair loss into a (K,) score."""
    X = Lm * w[None, :]
    if how == "mean":
        return X.mean(1)
    if how == "med":
        return np.median(X, axis=1)
    if how == "trim25":                      # drop the worst quarter of pairs
        Xs = np.sort(X, axis=1)
        return Xs[:, : max(1, int(0.75 * X.shape[1]))].mean(1)
    if how == "cvar25":                      # keep ONLY the worst quarter
        Xs = np.sort(X, axis=1)
        return Xs[:, max(1, int(0.75 * X.shape[1])):].mean(1)
    if how == "zset":                        # standardise each pair ACROSS the candidate set
        mu = Lm.mean(0); sg = np.maximum(Lm.std(0), 1e-9)
        return (((Lm - mu) / sg) * w[None, :]).mean(1)
    if how == "rankset":                     # nonparametric set-conditional
        r = np.argsort(np.argsort(Lm, axis=0, kind="stable"), axis=0, kind="stable")
        return ((r / max(Lm.shape[0] - 1, 1)) * w[None, :]).mean(1)
    raise KeyError(how)


# ------------------------------------------------------------------ the variant registry
#: (name, loss, weight, agg).  `shipped` must be first -- it is the reference column.
FACTORED = [
    ("shipped",        "risk", "ship",     "mean"),
    # --- the four shipped choices, ablated one at a time -------------------------
    ("risk_one",       "risk", "one",      "mean"),      # remove uncertainty weighting
    ("risk_ship2",     "risk", "ship2",    "mean"),      # sharpen it
    ("risk_shiphalf",  "risk", "shiphalf", "mean"),      # soften it
    ("risk_ent",       "risk", "ent",      "mean"),      # entropy instead of sd
    ("risk_maxp",      "risk", "maxp",     "mean"),
    ("risk_sep",       "risk", "sep",      "mean"),      # sequence-separation weighting
    ("risk_isep",      "risk", "isep",     "mean"),
    ("risk_lr",        "risk", "lr",       "mean"),      # long-range pairs only
    ("risk_sr",        "risk", "sr",       "mean"),      # local pairs only
    ("risk_shipsep",   "risk", "shipsep",  "mean"),
    ("risk_shipisep",  "risk", "shipisep", "mean"),
    ("risk_med",       "risk", "ship",     "med"),       # robust aggregations
    ("risk_trim25",    "risk", "ship",     "trim25"),
    ("risk_cvar25",    "risk", "ship",     "cvar25"),
    ("risk_zset",      "risk", "ship",     "zset"),      # set-conditional standardisation
    ("risk_rankset",   "risk", "ship",     "rankset"),
    ("risk_one_zset",  "risk", "one",      "zset"),
    # --- other residual forms ----------------------------------------------------
    ("abs_ship",       "abs",  "ship",     "mean"),
    ("abs_one",        "abs",  "one",      "mean"),
    ("z_one",          "z",    "one",      "mean"),      # z already carries 1/sd
    ("z_ship",         "z",    "ship",     "mean"),
    ("z_zset",         "z",    "one",      "zset"),
    ("z_med",          "z",    "one",      "med"),
    ("z_trim25",       "z",    "one",      "trim25"),
    ("z_cvar25",       "z",    "one",      "cvar25"),
    ("hub_one",        "hub",  "one",      "mean"),      # robust losses
    ("hub_ship",       "hub",  "ship",     "mean"),
    ("tuk_one",        "tuk",  "one",      "mean"),
    ("tuk_ship",       "tuk",  "ship",     "mean"),
    ("log_ship",       "log",  "ship",     "mean"),
    ("log_one",        "log",  "one",      "mean"),
    ("sq_ship",        "sq",   "ship",     "mean"),
    ("nll_one",        "nll",  "one",      "mean"),      # full likelihood, not Bayes risk
    ("nll_ship",       "nll",  "ship",     "mean"),
    ("nll_zset",       "nll",  "one",      "zset"),
    ("nll_lr",         "nll",  "lr",       "mean"),
]

BETAS = (0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.0)


def blend_scores(c, betas=BETAS, scale="pred", weight="ship", w=None):
    """The centre-blend family.  m(beta) = (1-beta)*expected + beta*setmean.

    beta = 0 is the distance objective in residual form; beta = 1 is TYPICALITY.
    `scale` selects the denominator: the predicted sd, or the set's own sd.
    """
    D = c["D"]
    s = c["sd"] if scale == "pred" else c["setsd"]
    ww = c["w"][weight] if w is None else w
    out = {}
    for b in betas:
        m = (1.0 - b) * c["exp"] + b * c["setmean"]
        out[f"blend{scale}_{b:.2f}"] = ((np.abs(D - m[None, :]) / np.maximum(s, 1e-3)[None, :])
                                        * ww[None, :]).mean(1)
    return out


def contact_scores(c):
    """Contact-vs-distance: collapse the prediction to P(d < t) and score the indicator."""
    centres, prob = c["prob"]
    out = {}
    for t in (8.0, 10.0):
        pc = prob[:, centres < t].sum(1)
        ind = (c["D"] < t).astype(float)
        # cross-entropy of the observed contact state under the predicted contact probability
        out[f"contact_ce{int(t)}"] = -(ind * np.log(np.maximum(pc, 1e-6))[None, :]
                                       + (1 - ind) * np.log(np.maximum(1 - pc, 1e-6))[None, :]).mean(1)
        out[f"contact_l1_{int(t)}"] = np.abs(ind - pc[None, :]).mean(1)
    return out


def spectral_scores(c, n):
    """Torgerson / low-rank form: match the eigenspectrum of the double-centred squared
    distance matrix.  A global shape comparison rather than a pairwise one."""
    K, P = c["D"].shape
    # rebuild full (n,n) squared-distance matrices; |i-j| < 2 entries are IDENTICAL in both
    # arms (ideal CA geometry) and so cannot contribute -- they are filled for both alike.
    def spec(Dv):
        M = np.zeros((len(Dv), n, n)) if Dv.ndim == 2 else np.zeros((1, n, n))
        Dv = np.atleast_2d(Dv)
        ii, jj = c["ii"], c["jj"]
        M[:, ii, jj] = Dv ** 2
        M[:, jj, ii] = Dv ** 2
        for d in (1,):
            k = np.arange(n - d)
            M[:, k, k + d] = 3.8 ** 2
            M[:, k + d, k] = 3.8 ** 2
        J = np.eye(n) - 1.0 / n
        G = -0.5 * np.einsum("ab,kbc,cd->kad", J, M, J)
        ev = np.linalg.eigvalsh(G)[:, ::-1][:, :3]
        return ev
    ep = spec(c["exp"])[0]
    ec = spec(c["D"])
    return {"spectral3": np.abs(ec - ep[None, :]).sum(1)}


def all_scores(p, c):
    s = {}
    for name, loss, w, agg in FACTORED:
        s[name] = _agg(c[loss], c["w"][w], agg)
    s.update(blend_scores(c, scale="pred"))
    s.update(blend_scores(c, scale="set"))
    s.update(contact_scores(c))
    c["ii"], c["jj"] = p["i"], p["j"]
    s.update(spectral_scores(c, p["n"]))
    return s


# ------------------------------------------------------------------ rho diagnostics
def _rank(x):
    return np.argsort(np.argsort(np.asarray(x, float), kind="stable"), kind="stable").astype(float)


def spearman(a, b):
    ra, rb = _rank(a), _rank(b)
    ra -= ra.mean(); rb -= rb.mean()
    den = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / den) if den > 0 else 0.0


def rhos(score, rr, band=1.5):
    """(within-target rho, in-band rho).  Positive means the score orders correctly."""
    rw = spearman(score, rr)
    m = rr <= rr.min() + band
    rb = spearman(score[m], rr[m]) if m.sum() >= 8 else float("nan")
    return rw, rb


# ------------------------------------------------------------------ the sweep
def run(K=500, targets=None, tag=None, verbose=True):
    tg = targets if targets is not None else L.targets()
    tag = tag or f"K{K}"
    rows = []
    t0 = time.time()
    for q, t in enumerate(tg):
        p = L.pack(t["pdb"], K)
        c = context(p)
        rr = p["rr"]
        sc = all_scores(p, c)
        rec = {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "k": int(p["k"]),
               "U": int(p["U"]), "best": float(rr.min()), "mean": float(rr.mean()),
               "median": float(np.median(rr)), "rand": L.rand_expect(rr), "sel": {}, "rho": {},
               "rho_band": {}, "ties": {}}
        for name, s in sc.items():
            rec["sel"][name] = L.sel_of(s, rr)
            rw, rb = rhos(s, rr)
            rec["rho"][name] = rw
            rec["rho_band"][name] = rb
            rec["ties"][name] = int((s <= s.min()).sum())
        rows.append(rec)
        if verbose and (q + 1) % 25 == 0:
            print(f"  {q+1}/{len(tg)}  ({time.time()-t0:.0f}s)", flush=True)
            json.dump({"K": K, "rows": rows}, open(os.path.join(L.RESULTS, f"sel_obj_{tag}.json"), "w"))
    json.dump({"K": K, "rows": rows}, open(os.path.join(L.RESULTS, f"sel_obj_{tag}.json"), "w"))
    return rows


# ------------------------------------------------------------------ reporting
def report(tag="K500"):
    d = json.load(open(os.path.join(L.RESULTS, f"sel_obj_{tag}.json")))
    rows = d["rows"]
    fold = np.array([r["fold"] for r in rows], int)
    names = list(rows[0]["sel"].keys())
    S = np.array([[r["sel"][v] for v in names] for r in rows])          # (T, V)
    RHO = np.array([[r["rho"][v] for v in names] for r in rows])
    RHB = np.array([[r["rho_band"][v] for v in names] for r in rows])
    best = np.array([r["best"] for r in rows])
    rand = np.array([r["rand"] for r in rows])
    ship = S[:, names.index("shipped")]
    n = len(rows)

    print(f"\n{'='*100}\nOBJECTIVE DECOMPOSITION  --  {tag},  n = {n} targets,  "
          f"identical candidate set per target\n{'='*100}")
    print(f"ORACLE CEILING (best member)     {best.mean():7.3f}")
    print(f"REALIZED, shipped objective      {ship.mean():7.3f}     GAP {ship.mean()-best.mean():+7.3f}")
    print(f"MATCHED RANDOM (exact expect.)   {rand.mean():7.3f}     shipped - random "
          f"{ship.mean()-rand.mean():+7.3f}")

    print(f"\nA. EVERY VARIANT, sorted by mean realised RMSD.  'vs shipped' is paired, "
          f"target as the unit.\n   rho_w = mean per-target Spearman (the design equation's "
          f"parameter); rho_b = in-band (<= best+1.5 A).")
    print(f"  {'variant':<22}{'RMSD':>8}{'median':>8}{'vs shipped [95% CI]':>28}{'W/L':>9}"
          f"{'rho_w':>8}{'rho_b':>8}")
    order = np.argsort(S.mean(0))
    for v in order:
        r = L.report_pair(names[v], S[:, v], ship, fold)
        print(f"  {names[v]:<22}{S[:,v].mean():>8.3f}{np.median(S[:,v]):>8.3f}"
              f"   {r['diff']:+.3f} [{r['lo']:+.3f},{r['hi']:+.3f}]{r['W']:>6}/{r['L']:<4}"
              f"{np.nanmean(RHO[:,v]):>8.3f}{np.nanmean(RHB[:,v]):>8.3f}")

    print(f"\nB. THE HEADLINE: LEAVE-FOLD-OUT VARIANT SELECTION  (PREREG E2)")
    lfo, choice, ins = L.lfo_variant(S, fold)
    print(f"  in-sample best-of-{len(names)}   {names[ins]:<22}{S[:,ins].mean():7.3f}   "
          f"[OVERFIT UPPER BOUND, not a result]")
    print(f"  leave-fold-out realised                        {lfo.mean():7.3f}")
    print(f"  per-fold choice: " + ", ".join(f"{k}:{names[v]}" for k, v in choice.items()))
    for nm, ref in (("LFO vs shipped", ship), ("LFO vs random", rand)):
        print("  " + L.fmt_pair(L.report_pair(nm, lfo, ref, fold)))
    print("  " + L.fmt_pair(L.report_pair("shipped vs random", ship, rand, fold)))

    print(f"\nC. THE CENTRE-BLEND FAMILY  (distance objective at beta=0 -> typicality at beta=1)")
    for scale in ("pred", "set"):
        print(f"  scale = {scale}")
        for b in BETAS:
            v = names.index(f"blend{scale}_{b:.2f}")
            r = L.report_pair("", S[:, v], ship, fold)
            print(f"    beta {b:4.2f}  RMSD {S[:,v].mean():7.3f}  vs shipped {r['diff']:+.3f} "
                  f"[{r['lo']:+.3f},{r['hi']:+.3f}]  {r['W']}W/{r['L']}L  "
                  f"rho_w {np.nanmean(RHO[:,v]):+.3f}  rho_b {np.nanmean(RHB[:,v]):+.3f}")

    print(f"\nD. WHAT THE SHIPPED WEIGHTING IS WORTH  (its own ablations)")
    for v in ("risk_one", "risk_ship2", "risk_shiphalf", "risk_ent", "risk_maxp"):
        r = L.report_pair(v, S[:, names.index(v)], ship, fold)
        print("  " + L.fmt_pair(r))
    return dict(names=names, S=S, RHO=RHO, RHB=RHB, best=best, rand=rand, fold=fold, rows=rows)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        report(sys.argv[2] if len(sys.argv) > 2 else "K500")
    else:
        K = int(sys.argv[1]) if len(sys.argv) > 1 else 500
        run(K=K, tag=f"K{K}")
        report(f"K{K}")
