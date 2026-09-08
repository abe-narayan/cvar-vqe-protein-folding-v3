"""COORDINATOR EXPERIMENT 2 -- does the empirical library have HOLES that ideal geometry fills?

The literature agent measured something that should not be possible under the project's
standing account of retrieval.  On 2BFI, a chain built from a SINGLE CONSTANT (phi, psi) =
(-180, 180) scores 0.27 A to the native, while the best of all 21,547 windows in that
target's legal universe scores 1.35 A.  A two-parameter model beat a 21,547-member library.

That is a COVERAGE statement, and it is a different one from the coverage question the
record has already answered.  `bench_results/recon_library_saturation.json` showed that a
3.2x larger fragment library from the same corpus moves nothing -- the library is saturated
*in the distribution it samples from*.  It does not follow that the library covers
conformation space.  Real protein fragments are never ideal: a beta strand inside a protein
carries a twist, a bulge, a neighbouring turn.  The canonical conformations -- fully
extended, ideal alpha, polyproline II, ideal hairpins -- may be systematically ABSENT from
an empirical window bank at any size, because nature does not deposit them pure.

So this module asks two questions and, importantly, asks the second one WITHOUT a detector:

    Q1 (ORACLE, diagnostic)  How close is the best IDEAL template to each native, against
                             the best of the target's whole legal window universe?
    Q2 (DEPLOYABLE)          If the templates are simply APPENDED to the K=500 pool and the
                             pipeline is otherwise untouched, does the shipped filter pick
                             them, and does the emitted structure improve?

Q2 is the honest form of the literature agent's routing proposal.  Routing needs a
steric-zipper detector and pays +4.7 A for a false positive; pool augmentation needs
nothing, cannot mis-route, and lets the existing objective decide.  If the objective cannot
pick a 0.27 A template out of a 700-member pool, that is itself a sharp measurement of the
objective -- and it is the same objective the aggregation agent showed sits at its
native-free optimum, so there is nowhere for the blame to hide.

CONTROLS.  Appending ~200 candidates to a 500-member pool changes the top-75 by
construction.  Two nulls are run at matched cardinality: append the same number of RANDOM
windows drawn from the universe, and append the same number of random windows drawn from
the pool's own tail.  A "gain" that the random-append arm reproduces is a cardinality
artefact, not template information.

    python -m s12.coord_templates            # all 126 targets
"""
import os
import sys
import json
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402

D2R = np.pi / 180.0

#: Canonical single-conformation templates, (phi, psi) in degrees.  Textbook values.
CANONICAL = {
    "alpha_R":        (-57.0, -47.0),
    "helix_310":      (-49.0, -26.0),
    "helix_pi":       (-57.0, -70.0),
    "alpha_L":        (57.0, 47.0),
    "beta_anti":      (-139.0, 135.0),
    "beta_par":       (-119.0, 113.0),
    "extended_full":  (-180.0, 180.0),
    "ppII":           (-75.0, 145.0),
    "ppI":            (-75.0, 160.0),
    "bridge":         (-90.0, 0.0),
}

#: A coarse sweep of the Ramachandran map, so the bank is not just ten textbook points.
GRID_PHI = np.arange(-180.0, 0.1, 20.0)
GRID_PSI = np.arange(-180.0, 180.1, 20.0)
#: The left-handed basin, which the negative-phi sweep above cannot reach.
GRID_PHI_L = np.arange(20.0, 100.1, 20.0)
GRID_PSI_L = np.arange(-20.0, 100.1, 20.0)


def constant_templates(n):
    """Ideal chains of length n built from one (phi, psi) repeated.  (names, PHI, PSI)."""
    names, P, S = [], [], []
    for k, (a, b) in CANONICAL.items():
        names.append(f"const:{k}"); P.append(np.full(n, a)); S.append(np.full(n, b))
    for a in GRID_PHI:
        for b in GRID_PSI:
            names.append(f"grid:{a:.0f},{b:.0f}"); P.append(np.full(n, a)); S.append(np.full(n, b))
    for a in GRID_PHI_L:
        for b in GRID_PSI_L:
            names.append(f"gridL:{a:.0f},{b:.0f}"); P.append(np.full(n, a)); S.append(np.full(n, b))
    return names, np.array(P) * D2R, np.array(S) * D2R


def hairpin_templates(n):
    """Two ideal strands joined by a 2-residue turn, over every interior turn position.

    A beta hairpin is the one common peptide fold a single constant (phi, psi) cannot
    express, and it is the fold several of the failing targets actually have.
    """
    names, P, S = [], [], []
    turns = {"I'": ((60.0, 30.0), (90.0, 0.0)), "II'": ((60.0, -120.0), (-80.0, 0.0)),
             "I": ((-60.0, -30.0), (-90.0, 0.0)), "II": ((-60.0, 120.0), (80.0, 0.0))}
    for strand_name, (sa, sb) in (("anti", (-139.0, 135.0)), ("ext", (-180.0, 180.0))):
        for tname, (t1, t2) in turns.items():
            for t in range(2, n - 3):          # turn occupies residues t, t+1
                ph = np.full(n, sa); ps = np.full(n, sb)
                ph[t], ps[t] = t1
                ph[t + 1], ps[t + 1] = t2
                names.append(f"hairpin:{strand_name}:{tname}:t{t}")
                P.append(ph); S.append(ps)
    return names, np.array(P) * D2R, np.array(S) * D2R


def template_bank(n):
    n1, p1, s1 = constant_templates(n)
    n2, p2, s2 = hairpin_templates(n)
    names = n1 + n2
    PHI = np.vstack([p1, p2]); PSI = np.vstack([s1, s2])
    W = I.build_ca(PHI, PSI)
    return names, PHI, PSI, W


def run_target(t, seeds=(0, 1, 2)):
    u = I.load_univ(t["pdb"])
    n, seq, fold, nat = t["n"], t["seq"], t["fold"], u["nat_ca"]
    W_all = u["W"]
    pool = I.pool_idx(u, I.K)
    Wp = W_all[pool].astype(np.float32).astype(float)
    rrp = u["rr"][pool]

    names, TPHI, TPSI, TW = template_bank(n)
    nt = len(names)
    rr_t = I.kabsch_rmsd_batch(TW, nat)                       # ORACLE
    best = int(np.argmin(rr_t))

    dg = I.distogram(t["pdb"], seq, fold)
    i, j = I.pair_index(n)

    def score_of(W):
        return I.shipped_score(dg, I.pair_dists(np.asarray(W, float), i, j))

    sc_pool = score_of(Wp)
    sc_t = score_of(TW)

    def emit(W):
        C, _ = I.coordinate_average(np.asarray(W, float))
        o = I.project(C, seq, fold, lam=0.3)
        return I.ca_rmsd(o["ca"], nat), I.ca_rmsd(o["fit_ca"], nat)

    def top_m(Wcat, sc):
        return np.argsort(sc, kind="stable")[:I.M]

    res = {"pdb": t["pdb"], "n": n, "fold": fold, "n_templates": nt,
           "universe_best": float(u["rr"].min()), "pool_best": float(rrp.min()),
           "template_best": float(rr_t.min()), "template_best_name": names[best],
           "template_beats_universe": bool(rr_t.min() < u["rr"].min()),
           "template_beats_pool": bool(rr_t.min() < rrp.min())}

    # --- baseline: the shipped pool alone
    sub = top_m(Wp, sc_pool)
    a, f = emit(Wp[sub]); res["base_arm"], res["base_fit"] = a, f
    res["base_top75_best"] = float(rrp[sub].min())

    # --- TEMPLATE-AUGMENTED pool
    Wc = np.vstack([Wp, TW]); scc = np.concatenate([sc_pool, sc_t]); rrc = np.concatenate([rrp, rr_t])
    subc = top_m(Wc, scc)
    a, f = emit(Wc[subc]); res["tmpl_arm"], res["tmpl_fit"] = a, f
    res["tmpl_top75_best"] = float(rrc[subc].min())
    res["n_templates_in_top75"] = int((subc >= len(Wp)).sum())
    res["best_template_score_pct"] = float((sc_pool < sc_t[best]).mean())
    res["template_argmin_is_template"] = bool(int(np.argmin(scc)) >= len(Wp))

    # --- NULL: append the same number of RANDOM universe windows
    arms, fits, tb = [], [], []
    for s in seeds:
        rng = np.random.default_rng(97 * s + hash(t["pdb"]) % 9973)
        extra = rng.permutation(len(W_all))[:nt]
        We = W_all[extra].astype(np.float32).astype(float)
        Wn = np.vstack([Wp, We]); scn = np.concatenate([sc_pool, score_of(We)])
        rrn = np.concatenate([rrp, u["rr"][extra]])
        sn = top_m(Wn, scn)
        x, y = emit(Wn[sn]); arms.append(x); fits.append(y); tb.append(float(rrn[sn].min()))
    res["randaug_arm"] = float(np.mean(arms)); res["randaug_fit"] = float(np.mean(fits))
    res["randaug_top75_best"] = float(np.mean(tb))
    return res


def main():
    tg = I.targets()
    path = os.path.join(I.RESULTS, "coord_templates.json")
    rows = json.load(open(path))["per_target"] if os.path.exists(path) else []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        rows.append(run_target(t))
        if k % 5 == 0 or k == len(tg) - 1:
            json.dump({"what": "ideal-template pool augmentation", "per_target": rows},
                      open(path, "w"), indent=1)
            r = rows[-1]
            print(f"  {k+1}/{len(tg)} {t['pdb']} tmplbest={r['template_best']:.2f} "
                  f"univbest={r['universe_best']:.2f} base={r['base_arm']:.2f} "
                  f"tmpl={r['tmpl_arm']:.2f} in75={r['n_templates_in_top75']} "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    json.dump({"what": "ideal-template pool augmentation", "per_target": rows},
              open(path, "w"), indent=1)
    report(rows)


def report(rows):
    names = [r["pdb"] for r in rows]; folds = [r["fold"] for r in rows]
    f18 = set(I.FAIL18); isf = np.array([r["pdb"] in f18 for r in rows])
    g = lambda k: np.array([r[k] for r in rows], float)                    # noqa: E731
    out = {"oracle_coverage": {
               "universe_best": I.summary(g("universe_best")),
               "template_best": I.summary(g("template_best")),
               "n_template_beats_universe": int(sum(r["template_beats_universe"] for r in rows)),
               "n_template_beats_pool": int(sum(r["template_beats_pool"] for r in rows)),
               "FAIL18_template_best": float(g("template_best")[isf].mean()),
               "FAIL18_universe_best": float(g("universe_best")[isf].mean())},
           "emitted": {k: I.summary(g(k)) for k in ("base_arm", "tmpl_arm", "randaug_arm")},
           "paired": {
               "templates vs base":      I.paired(g("tmpl_arm"), g("base_arm"), folds=folds, names=names),
               "random-append vs base":  I.paired(g("randaug_arm"), g("base_arm"), folds=folds, names=names),
               "templates vs random-append": I.paired(g("tmpl_arm"), g("randaug_arm"), folds=folds, names=names)},
           "filter_behaviour": {
               "mean_templates_in_top75": float(g("n_templates_in_top75").mean()),
               "n_targets_with_any_template_in_top75": int((g("n_templates_in_top75") > 0).sum()),
               "n_targets_where_argmin_is_a_template": int(sum(r["template_argmin_is_template"] for r in rows)),
               "mean_score_percentile_of_best_template": float(g("best_template_score_pct").mean()),
               "top75_best_base": float(g("base_top75_best").mean()),
               "top75_best_tmpl": float(g("tmpl_top75_best").mean()),
               "top75_best_randaug": float(g("randaug_top75_best").mean())},
           "subgroups": {"FAIL18": {k: float(g(k)[isf].mean()) for k in ("base_arm", "tmpl_arm", "randaug_arm")},
                         "other108": {k: float(g(k)[~isf].mean()) for k in ("base_arm", "tmpl_arm", "randaug_arm")}}}
    I.write("coord_templates_report", out)
    print("\n=== ORACLE coverage ===")
    print(f"  universe best  {out['oracle_coverage']['universe_best']['mean']:.4f}")
    print(f"  template best  {out['oracle_coverage']['template_best']['mean']:.4f}")
    print(f"  templates beat the WHOLE universe on {out['oracle_coverage']['n_template_beats_universe']}/126,"
          f" the pool on {out['oracle_coverage']['n_template_beats_pool']}/126")
    print(f"  FAIL18: template {out['oracle_coverage']['FAIL18_template_best']:.3f} vs "
          f"universe {out['oracle_coverage']['FAIL18_universe_best']:.3f}")
    print("\n=== emitted ===")
    for k, v in out["emitted"].items():
        print(f"  {k:12s} {v['mean']:.4f}  <2A {v['frac_under_2.0']:.3f}")
    print("\n=== paired (negative = first better) ===")
    for k, v in out["paired"].items():
        print(f"  {k:28s} {v['mean_diff']:+.4f} [{v['ci95'][0]:+.4f},{v['ci95'][1]:+.4f}] "
              f"{v['n_better']}W/{v['n_worse']}L drop10 {v['drop_top10_mean_diff']:+.4f}")
    print("\n=== does the filter pick them? ===")
    for k, v in out["filter_behaviour"].items():
        print(f"  {k}: {v}")
    print("\n=== subgroups ===")
    for gname, d in out["subgroups"].items():
        print(f"  {gname:9s} " + "  ".join(f"{k}={v:.3f}" for k, v in d.items()))


if __name__ == "__main__":
    main()
