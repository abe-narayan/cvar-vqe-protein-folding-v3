"""s16/integrate.py -- THE INTEGRATED ARCHITECTURE, and the Legacy/AMBER ablation INSIDE it.

    CVaR-VQE ensemble  ->  structure-level readout  ->  Legacy detection  ->  AMBER repair

This is the one experiment in which all four mandated pillars meet on the same targets, the
same ensembles and the same readout, so that Legacy and AMBER are EXPLICITLY COMPARABLE: they
are applied as *the same operation* -- drop the worst quartile of the ensemble by that energy
-- differing only in which energy does the ranking.  Anything else compares an energy to an
operation.

WHAT IS GENUINE HERE.  The VQE is `s14.vqe_run.run`: a real parameterised circuit (MPS ansatz)
sampled by shots, optimised by Adam against a CVaR tail of a real qubit Hamiltonian, under a
hard objective-evaluation budget (`s14.vqe_lib.Counter`).  `alpha = 1.0` is plain VQE and
`alpha < 1` is CVaR-VQE, so both pillars are present as arms of one ladder rather than as
separate claims.  Legacy is the genuine eleven-component Miyazawa-Jennings-style model at
`DEFAULT_WEIGHTS`, never fitted.  AMBER is a genuine ff14SB/GBn2 single point computed on
demand (`s13.qarch_lib.amber_energies`, ~6 ms per configuration) -- NOT the 1.1% precomputed
table, which a sampled ensemble would almost never intersect.

THE DESIGN.

  generators (matched on objective evaluations, the only budget convention used here)
      cvar{a}   CVaR-VQE at alpha in {0.1, 0.25, 0.5, 1.0};  1.0 is plain VQE
      bestofN   the UNTRAINED circuit sampled to the same budget -- the mandatory control.
                The programme's standing result is that running the VQE is WORSE than not
                running it (0/12), and that result was only visible against THIS control.
      uniform   uniform random configurations, same budget

  ensemble  the m distinct sampled configurations of lowest objective, m in {5, 20, 75}

  readouts, all on the SAME ensemble, differing only in the drop rule
      ctrl      coordinate average of the members
      leg       drop the worst ceil(m/4) by genuine Legacy total, then average
      amb       drop the worst ceil(m/4) by genuine AMBER single point, then average
      legamb    Legacy then AMBER, dropping ceil(m/4) IN TOTAL so the count matches `leg`
                and `amb` and the interaction is not confounded with ensemble size
      legamb_f  Legacy then AMBER, each at full quartile strength (keeps ~9/16) -- the
                architecture as specified, reported but never used for INTER
      rnd       drop a uniformly random ceil(m/4) (2 draws) -- a filter must beat this
      tors      the SAME ensemble averaged in TORSION space instead

      LEGACY = leg - ctrl     AMBER = amb - ctrl     INTER = (legamb - amb) - (leg - ctrl)

WHY `tors` IS IN THE TABLE.  It is not a spare arm.  `s16/steer.py` established that RMSD is
strongly convex in torsion displacement -- half-way to the exact native torsions is worth 5%
of the answer -- so this architecture's choice to read out in coordinate space is a claim, and
`tors` is the measurement that prices it on this instrument.

HONESTY ABOUT n.  The enumerated instrument is NINE targets.  That is the price of a certified
optimum and an exhaustively enumerable space, and it is far below the 126-target instrument
every other number in this sprint is quoted on.  Nine is also inside the range where this
programme has had conclusions reverse -- most recently an n = 8 arm reading -0.106 with a CI
EXCLUDING zero where the n = 126 truth was -0.003 [-0.011, +0.007].  Every table here is
therefore printed PER TARGET as well as aggregated, win/loss counts are given, and no
aggregate is quoted without them.

ORACLE.  `Enum.rmsd` and the native trace are read only to score.  No native quantity enters a
generator, an ensemble, a drop rule, or a hyperparameter.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I              # noqa: E402
from s14 import vqe_lib as V                 # noqa: E402
from s14 import vqe_run as R                 # noqa: E402
from s15 import seed as SD                   # noqa: E402
from s16 import energy_lib as EL             # noqa: E402

ALPHAS = (0.1, 0.25, 0.5, 1.0)
MS = (5, 20, 75)
BUDGET = 6144
SHOTS = 512
SEEDS = (0, 1, 2)


def _ca_of(en, idx):
    """(B,) configuration indices -> (B, n, 3) ideal-geometry CA traces."""
    S = en.states(idx)
    rows = np.arange(en.n)
    phi = en.PHI[rows[None, :], S]
    psi = en.PSI[rows[None, :], S]
    return I.build_ca(phi, psi), phi, psi


def _untrained(en, budget, shots, seed):
    """best-of-N from the UNTRAINED circuit at the same budget. The mandatory control."""
    rng = np.random.default_rng(seed + 4441)
    an = R.make_ansatz("mps2f", en.n_qubits)
    th0 = R.init_theta(an, rng, 0.8, "random", None, bits_per_res=en.bits_per_res)
    seen = []
    for _ in range(max(1, budget // shots)):
        b = an.sample(th0, shots, rng)
        seen.append(R._bits_to_index(b))
    return np.concatenate(seen)[:budget]


def _uniform(en, budget, seed):
    return SD.stable_rng(en.pdb, f"unif{seed}").integers(0, en.N, budget)


def _ensembles(en, seed):
    """every generator's sampled index multiset, matched on objective evaluations."""
    out = {}
    for a in ALPHAS:
        r = R.run(en.prior, en.n_qubits, a, BUDGET, shots=SHOTS, seed=seed,
                  bits_per_res=en.bits_per_res, exact_dist=False, keep_seen=True)
        out[f"cvar{a}"] = np.asarray(r["seen"], np.int64)
    out["bestofN"] = _untrained(en, BUDGET, SHOTS, seed)
    out["uniform"] = _uniform(en, BUDGET, seed)
    return out


def _pick(en, seen, m):
    """the m DISTINCT sampled configurations of lowest objective."""
    u = np.unique(np.asarray(seen, np.int64))
    if u.size <= m:
        return u
    return u[np.argsort(en.prior[u], kind="stable")[:m]]


def _readouts(en, idx, nat, amber, legacy, rng):
    """every drop rule on the SAME ensemble. Returns {name: rmsd}."""
    W, phi, psi = _ca_of(en, idx)
    m = len(idx)
    drop = max(1, math.ceil(m / 4.0))
    keep = m - drop

    def avg(sel):
        if len(sel) == 0:
            return float("nan")
        if len(sel) == 1:
            return float(I.ca_rmsd(W[sel][0], nat))
        c, _b = I.coordinate_average(W[sel])
        return float(I.ca_rmsd(c, nat))

    allsel = np.arange(m)
    o = {"ctrl": avg(allsel)}
    o["leg"] = avg(np.argsort(legacy, kind="stable")[:keep])
    o["amb"] = avg(np.argsort(amber, kind="stable")[:keep])
    #: Legacy first, then AMBER among the survivors -- the ORDERED composition the
    #: architecture specifies (detect with Legacy, repair/refine with AMBER).  TWO versions,
    #: because they answer different questions and conflating them would corrupt the
    #: factorial:
    #:   `legamb`      drops the SAME ceil(m/4) in total, half by each energy, so all three
    #:                 filtered cells keep the same count and the INTER contrast is a clean
    #:                 interaction rather than a confound with ensemble size;
    #:   `legamb_full` applies each filter at FULL quartile strength, keeping ~9/16 -- the
    #:                 architecture as actually specified, but NOT count-matched, so it is
    #:                 reported beside the others and never used for INTER.
    d1 = max(1, drop // 2)
    lk = np.argsort(legacy, kind="stable")[:m - d1]
    o["legamb"] = avg(lk[np.argsort(amber[lk], kind="stable")[:keep]])
    leg_keep = np.argsort(legacy, kind="stable")[:keep]
    k2 = max(1, len(leg_keep) - max(1, math.ceil(len(leg_keep) / 4.0)))
    o["legamb_full"] = avg(leg_keep[np.argsort(amber[leg_keep], kind="stable")[:k2]])
    rr = []
    for _ in range(2):
        rr.append(avg(rng.permutation(m)[:keep]))
    o["rnd"] = float(np.mean(rr))
    #: the torsion-space readout, priced rather than assumed away (see the module docstring).
    o["tors"] = float(I.ca_rmsd(
        I.build_ca(np.arctan2(np.sin(phi).mean(0), np.cos(phi).mean(0)),
                   np.arctan2(np.sin(psi).mean(0), np.cos(psi).mean(0))), nat))
    #: two reference points that are not readouts: what the ensemble contained.
    o["_best_member"] = float(en.rmsd[idx].min())          # ORACLE
    o["_mean_member"] = float(en.rmsd[idx].mean())         # ORACLE
    return o


def run(targets=None, seeds=SEEDS):
    from s13 import qarch_lib as QL
    tg = list(targets or V.ENUM_TARGETS)
    rows = []
    t0 = time.time()
    for ti, pdb in enumerate(tg):
        en = V.Enum(pdb)
        sp = QL.Space(pdb, en.k, seq=en.seq, n=en.n, fold=en.fold)
        nat = np.asarray(sp.nat_ca, float)
        for seed in seeds:
            gens = _ensembles(en, seed)
            for gname, seen in gens.items():
                for m in MS:
                    idx = _pick(en, seen, m)
                    S = en.states(idx)
                    #: GENUINE ff14SB/GBn2 single points, computed on demand.  The
                    #: precomputed table covers 1.1% of the space and a sampled ensemble
                    #: would essentially never intersect it.
                    amber = np.asarray(QL.amber_energies(sp, S), float)
                    comp = EL.legacy_components_of_windows(
                        en.seq, en.PHI[np.arange(en.n)[None, :], S],
                        en.PSI[np.arange(en.n)[None, :], S])
                    legacy = np.asarray(EL.legacy_total_from(comp), float)
                    rng = SD.stable_rng(pdb, f"int{seed}{gname}{m}")
                    o = _readouts(en, idx, nat, amber, legacy, rng)
                    rows.append({"pdb": pdb, "seed": int(seed), "gen": gname,
                                 "m": int(m), "n_ens": int(len(idx)), **o})
        json.dump({"rows": rows}, open(os.path.join(RESULTS, "integrate.json"), "w"))
        print(f"  {ti+1}/{len(tg)} {pdb} done  ({time.time()-t0:.0f}s)", flush=True)
    report(rows)
    return rows


def _boot(d, rng, B=4000):
    d = np.asarray(d, float); k = len(d)
    if k == 0:
        return float("nan"), float("nan"), float("nan")
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def report(rows):
    rng = SD.stable_rng("integrate", "report")
    gens = ["cvar0.1", "cvar0.25", "cvar0.5", "cvar1.0", "bestofN", "uniform"]
    pdbs = sorted({r["pdb"] for r in rows})
    nseed = len({r["seed"] for r in rows})
    print(f"\nn = {len(pdbs)} ENUMERATED targets, {nseed} seeds. "
          f"This is NOT the 126-target instrument; read the per-target table too.\n")

    for m in MS:
        print(f"=== ensemble size m = {m} " + "=" * 46)
        print(f"  {'generator':<12}{'ctrl':>8}{'leg':>8}{'amb':>8}{'legamb':>8}"
              f"{'lg+ab_f':>8}{'rnd':>8}{'tors':>8}   {'LEGACY - rnd [95% CI]':>26}")
        for g in gens:
            sel = [r for r in rows if r["gen"] == g and r["m"] == m]
            if not sel:
                continue
            mu = {k: np.mean([r[k] for r in sel]) for k in
                  ("ctrl", "leg", "amb", "legamb", "legamb_full", "rnd", "tors")}
            d = [r["leg"] - r["rnd"] for r in sel]
            a, lo, hi = _boot(d, rng)
            print(f"  {g:<12}" + "".join(f"{mu[k]:>8.3f}" for k in
                  ("ctrl", "leg", "amb", "legamb", "legamb_full", "rnd", "tors"))
                  + f"   {a:+.3f} [{lo:+.3f},{hi:+.3f}]")
        print()

    print("=== the factorial effects, pooled over generators, per ensemble size " + "=" * 8)
    print(f"  {'m':>4}{'LEGACY = leg-ctrl':>26}{'AMBER = amb-ctrl':>26}{'INTER':>26}")
    for m in MS:
        sel = [r for r in rows if r["m"] == m]
        L, Lo, Lh = _boot([r["leg"] - r["ctrl"] for r in sel], rng)
        A, Ao, Ah = _boot([r["amb"] - r["ctrl"] for r in sel], rng)
        X, Xo, Xh = _boot([(r["legamb"] - r["amb"]) - (r["leg"] - r["ctrl"]) for r in sel], rng)
        print(f"  {m:>4}{L:>+14.3f} [{Lo:+.3f},{Lh:+.3f}]{A:>+14.3f} [{Ao:+.3f},{Ah:+.3f}]"
              f"{X:>+14.3f} [{Xo:+.3f},{Xh:+.3f}]")

    print("\n=== coordinate vs torsion readout, the architecture's own design claim " + "=" * 5)
    for m in MS:
        sel = [r for r in rows if r["m"] == m]
        d = [r["tors"] - r["ctrl"] for r in sel]
        a, lo, hi = _boot(d, rng)
        w = sum(1 for x in d if x > 0); l = sum(1 for x in d if x < 0)
        print(f"  m = {m:<4} torsion minus coordinate  {a:+.3f} [{lo:+.3f},{hi:+.3f}]"
              f"   {w}W/{l}L for coordinate")

    print("\n=== PER TARGET, m = 75, pooled over seeds  (nine targets is not a population) ===")
    print(f"  {'pdb':<8}" + "".join(f"{g:>11}" for g in gens))
    for p in pdbs:
        line = f"  {p:<8}"
        for g in gens:
            s = [r["ctrl"] for r in rows if r["pdb"] == p and r["gen"] == g and r["m"] == 75]
            line += f"{np.mean(s):>11.3f}" if s else f"{'-':>11}"
        print(line)


if __name__ == "__main__":
    run()
