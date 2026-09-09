"""S24 / LANE D -- THE SET-EQUALITY DETERMINATION, SETTLED FROM SOURCE AND CONFIRMED NUMERICALLY.

THIS IS NOT A DIRECTIONAL EXPERIMENT. It reads no RMSD, touches no pool, and has no
outcome that can favour a hypothesis: it is a PROPERTY TEST of two library functions,
executed on adversarial synthetic inputs. Rule 0's six-axis fork enumeration therefore
does not apply and is not claimed -- what IS claimed is stated in `s24/PREREG_D.md` as a
falsifiable property with an explicit falsifier ("exhibit ONE (E, p, alpha) triple whose
CVaR tail support is not an energy-ordered prefix").

THE QUESTION (BRIEF Lane D): is the Sprint-22 set-equality result a theorem about the
ALGORITHM -- in which case it survives any candidate manifold and Sprint 24 must not
spend compute trying to break it -- or an artefact of THIS pool's energy structure?

THE ANSWER, from `core/quantum.py`:  IT IS A THEOREM ABOUT THE ALGORITHM.

  (1) `tail_indices(energies, alpha)`  [quantum.py:229-272]
      k = ceil(alpha * n)                                   -- line 252, depends on (n, alpha) ONLY
      order = np.argsort(e, kind="stable"); tail[order[:k]]  -- lines 256-259
      part = np.argpartition(e, k-1); tail = e < q; + ties   -- lines 261-271
      The signature takes NO probability vector and NO parameters. The returned mask is a
      pure function of the ENERGY ARRAY. The tail is the k lowest energies, by construction,
      for every input. There is no pool of any composition for which this is false.

  (2) `cvar_from_probs(energies, probs, alpha)`  [quantum.py:336-359] -- the one place
      p_theta actually enters the tail:
          order = np.argsort(e, kind="stable")               -- line 352  ENERGY order
          cum   = np.cumsum(p[order])                        -- line 353
          take  = np.clip(alpha - (cum - p[order]), 0, p[order])  -- line 354
      `take[j] > 0` requires BOTH conjuncts of that clip:
          (i)  `alpha - (cum[j] - p[order][j]) > 0`, i.e. the mass STRICTLY BEFORE position j
               in the ENERGY order is < alpha.  `cum - p[order]` is the exclusive prefix sum,
               non-decreasing in j because p >= 0, so it crosses alpha exactly once and (i)
               alone defines an initial PREFIX of the energy order, ending at a cut `m_cut`.
          (ii) `p[order][j] > 0`.  A zero-probability state inside that prefix is PUNCHED OUT.

CORRECTION OF RECORD, 2026-09-08 (coordinator, and it is right; the wording below is
the coordinator's verbatim and supersedes this file's first draft):

  > The realised CVaR tail's support is always a SUBSET of an initial prefix of the energy
  > order, and equals that prefix exactly when every state in the prefix carries positive
  > probability. p_theta can delete a member; it can never add one outside the classical
  > top-m.

  My first draft claimed the support is always an INITIAL prefix, dropping conjunct (ii).
  That is false: with exact zeros the support is a prefix WITH HOLES. My original 2916-cell
  test reported zero failures only because every probability family in it carried an epsilon
  floor (`+1e-30`, `1e-18`, `+1e-6`), so no exact zero was ever produced and the assertion
  could not fire. Both the assertion and the families are fixed below; the failure mode is
  now reproduced deliberately (`exact_zeros_*`) so the test has a live falsifier.

  EQUALITY IS THE EMPIRICAL REGIME; SUBSET-HOOD IS THE THEOREM. A trained RY/CNOT state has
  full support at generic angles, which is exactly why s22's Gate 1 passed at 100% on every
  cell that was actually run: with p > 0 everywhere, conjunct (ii) never bites and subset-hood
  collapses to equality.

  CONSEQUENCE, and the load-bearing half is UNCHANGED and if anything stronger. The trained
  state can move exactly TWO things and no third:
      (a) WHERE the prefix cuts, i.e. m -- a position on the classical top-m ladder -- and,
          where it has exact zeros, WHICH members of that prefix it DELETES;
      (b) the WEIGHTS inside the prefix -- closed by s23 L8 at four temperatures.
  It can NEVER ADD a candidate the classical energy order excluded. Deletion is a coverage
  loss with no information channel attached: the deleted member is chosen by the amplitude
  pattern, which is a function of the SAME energies the classical rank order already uses.
  So membership is bounded above by `argsort(E)[:m]` as an IDENTITY, not as an empirical
  finding -- independent of the landscape's shape, degeneracy, multimodality and scale, and
  therefore of the candidate manifold that produced it. CHANGING THE POOL CANNOT BREAK IT.

SCOPE -- the two (and only two) conditions under which set equality can fail, both of
which are degradations rather than channels:

  S1. DIFFERENT ENERGY. If the VQE's H is not the functional the classical control ranks
      with, the two sets differ trivially. That is a different-Hamiltonian comparison, not
      a quantum effect, and the harness holds H fixed across both arms so it cannot occur
      by accident.
  S2. INCOMPLETE SUPPORT -- and this is the SAME FACT as conjunct (ii) above, arriving by a
      second route rather than a separate scope condition. Exactly-zero probability deletes a
      prefix member on the EXACT path; failing to draw a candidate deletes it on the SAMPLED
      path (`ansatz.sample` -> `cvar` / `tail_indices`, where the tail is the bottom-k of the
      SAMPLED multiset). Either way the realised tail is a SUBSET of the classical prefix and
      p_theta can only DELETE, never ADD. One statement, two mechanisms: coverage loss with no
      information channel. Any deviation on a sampled path is a shot-budget artefact to be
      reported as such, not a result.

WHAT THIS FILE ADDS to the source reading: an adversarial numerical confirmation, so the
claim is measured and not merely argued. Nine energy structures chosen to be exactly the
things a "different candidate manifold" could plausibly change -- heavy degeneracy, exact
ties, bimodality, extreme scale, near-constant, inverted, discrete, one outlier, and a
pathological all-equal case -- crossed with nine probability structures including ones
deliberately concentrated on the HIGHEST-energy states (the adversarial case: if a pool's
energy structure could break the theorem, an adversarial p is where it would show).
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import quantum as Q          # noqa: E402
from s15 import seed as SD             # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)


# ------------------------------------------------------------------ energy structures
def energy_families(n: int, rng) -> dict:
    """Nine energy landscapes. Each is a thing a changed candidate manifold could produce."""
    return {
        "continuous_gaussian": rng.normal(0, 1, n),
        "heavy_degenerate": np.round(rng.normal(0, 1, n), 1),          # many exact ties
        "all_tied": np.zeros(n),                                       # total degeneracy
        "bimodal": np.where(rng.random(n) < 0.5, rng.normal(-5, .3, n),
                            rng.normal(+5, .3, n)),
        "huge_scale": rng.normal(0, 1, n) * 1e7,
        "near_constant": 1.0 + rng.normal(0, 1e-12, n),
        "one_outlier": np.r_[rng.normal(0, 1, n - 1), 1e9],
        "discrete_few": rng.integers(0, 3, n).astype(float),
        "monotone_ramp": np.arange(n, dtype=float),
    }


def prob_families(n: int, rng) -> dict:
    """Twelve probability vectors, adversarial in both directions.

    THE FIRST DRAFT'S BUG, kept visible: every family below that is meant to be sparse used
    to carry an epsilon floor (`+1e-30`, `1e-18`), so `p > 0` held everywhere and conjunct
    (ii) of `take`'s clip could never bite. The `exact_zeros_*` families use TRUE zeros and
    are what make the prefix-hood falsifier live. `float_floor_*` keep the old
    epsilon-floored versions so the difference between the two regimes is measured, not
    argued -- full support is the empirical regime a trained circuit is actually in.
    """
    idx = np.arange(n)
    u = rng.random(n) + 1e-6
    return {
        # --- FULL SUPPORT: the regime a trained RY/CNOT state is actually in
        "uniform": np.ones(n),
        "random": u,
        "sharp_random": rng.random(n) ** 12 + 1e-15,
        "exp_decay": np.exp(-idx / 3.0),
        "exp_growth": np.exp(+idx / 3.0),
        "float_floor_point_first": np.r_[1.0, np.full(n - 1, 1e-18)],
        "float_floor_point_last": np.r_[np.full(n - 1, 1e-18), 1.0],
        # --- EXACT ZEROS: holes are possible here, and this is where the old test failed
        "exact_zeros_random10pct": np.where(rng.random(n) < 0.10, 1.0, 0.0) + (idx == 0),
        "exact_zeros_random50pct": np.where(rng.random(n) < 0.50, 1.0, 0.0) + (idx == 0),
        "exact_zeros_alternating": (idx % 2 == 0).astype(float),
        "exact_zeros_point_last": np.where(idx == n - 1, 1.0, 0.0),
        "exact_zeros_two_spikes": np.isin(idx, [1, n - 2]).astype(float),
    }


def prefix_cut(e: np.ndarray, support: np.ndarray) -> np.ndarray:
    """The initial energy-order prefix that `support` must live inside.

    The prefix runs to the LAST supported position in the energy order. Everything the
    theorem asserts is about this set: `support` never escapes it, and the states inside it
    that are missing from `support` are exactly the zero-probability ones.
    """
    order = np.argsort(e, kind="stable")
    pos = np.flatnonzero(support[order])
    m_cut = int(pos[-1]) + 1 if pos.size else 0
    pref = np.zeros_like(support)
    pref[order[:m_cut]] = True
    return pref


def is_energy_ordered_prefix(e: np.ndarray, support: np.ndarray) -> bool:
    """PREFIX-HOOD (the strong, and FALSE-in-general, claim): support is an initial segment.

    Ties make the order non-unique, so the test is value-based rather than index-based:
    every supported state's energy must be <= every unsupported state's energy. This is the
    assertion the first draft made and it DOES fire on the `exact_zeros_*` families -- it is
    kept so the failure rate is a measured number rather than a concession.
    """
    if support.sum() == 0 or support.all():
        return True
    return float(e[support].max()) <= float(e[~support].min()) + 1e-12


def main() -> int:
    rng = SD.stable_rng("s24", "laneD", "setequality", 0, salt="s24d")
    rows, fails = [], []
    n_prefix_violations = n_full_support = n_full_support_equal = 0
    for n in (8, 16, 63, 64, 100, 512):
        E = energy_families(n, rng)
        P = prob_families(n, rng)
        for ename, e in E.items():
            for pname, p in P.items():
                for alpha in (0.02, 0.05, 0.15, 0.3, 0.5, 1.0):
                    # --- (2) the probability-weighted path, the ONLY one p_theta enters
                    val, q, mass = Q.cvar_from_probs(e, p, alpha)
                    sup = mass > 0
                    pref = prefix_cut(e, sup)

                    # === THE THEOREM (asserted): support NEVER escapes the energy prefix.
                    # p_theta can delete a member; it can never add one the classical order
                    # excluded. This is the load-bearing claim and it must never fail.
                    ok_subset = bool((sup & ~pref).sum() == 0)
                    # === and every hole in that prefix is EXACTLY a zero-probability state,
                    # which is conjunct (ii) of the clip, checked directly.
                    holes = pref & ~sup
                    ok_holes_are_zeros = bool(np.all(np.asarray(p, float)[holes] == 0.0)) \
                        if holes.any() else True

                    # === PREFIX-HOOD (the first draft's over-claim): measured, not asserted.
                    ok_prefix = is_energy_ordered_prefix(e, sup)
                    if not ok_prefix:
                        n_prefix_violations += 1
                    # === the FULL-SUPPORT corollary: where p > 0 everywhere -- the regime a
                    # trained circuit is actually in -- subset-hood collapses to EQUALITY.
                    full_support = bool(np.all(np.asarray(p, float) > 0.0))
                    eq_full = None
                    if full_support:
                        n_full_support += 1
                        m = int(sup.sum())
                        top = np.zeros(n, bool)
                        top[np.argsort(e, kind="stable")[:m]] = True
                        eq_full = bool(np.sort(e[sup]).tolist() == np.sort(e[top]).tolist())
                        n_full_support_equal += int(eq_full)

                    # --- (1) the index path: tail is a pure function of e, no p at all
                    k, qq, tail = Q.tail_indices(e, alpha)
                    ok_tail = is_energy_ordered_prefix(e, tail)
                    ok_size = (int(tail.sum()) == max(1, int(math.ceil(alpha * n))))
                    _, _, tail2 = Q.tail_indices(e, alpha)
                    ok_indep = bool((tail == tail2).all())

                    row = dict(n=n, energy=ename, probs=pname, alpha=alpha,
                               m=int(sup.sum()), m_cut=int(pref.sum()),
                               n_holes=int(holes.sum()), full_support=full_support,
                               ok_subset=ok_subset, ok_holes_are_zeros=ok_holes_are_zeros,
                               ok_prefix=bool(ok_prefix), eq_under_full_support=eq_full,
                               ok_tail=bool(ok_tail), ok_size=bool(ok_size),
                               ok_indep=ok_indep)
                    rows.append(row)
                    # ONLY the theorem-level assertions are failures. Prefix-hood is reported.
                    if not (ok_subset and ok_holes_are_zeros and ok_tail and ok_size
                            and ok_indep and (eq_full is not False)):
                        fails.append(row)

    # ---- the sharpest single statement, now run in the FULL-SUPPORT regime (a trained
    # ---- RY/CNOT state has generic angles and therefore p > 0 everywhere): vary p over
    # ---- 4000 random draws at FIXED e and count DISTINCT tail sets per realised size.
    n = 128
    e = rng.normal(0, 1, n)
    seen = {}
    for _ in range(4000):
        p = rng.random(n) ** rng.integers(1, 20)
        p = p / p.sum()
        assert (p > 0).all()
        _, _, mass = Q.cvar_from_probs(e, p, 0.15)
        sup = tuple(np.flatnonzero(mass > 0).tolist())
        seen.setdefault(len(sup), set()).add(sup)
    distinct_by_m = {int(k): len(v) for k, v in sorted(seen.items())}
    multi = {k: v for k, v in distinct_by_m.items() if v != 1}

    ok = (not fails) and (not multi)
    out = dict(
        experiment="s24_laneD_set_equality_property_test",
        revision=2,
        revision_note=(
            "REV2 2026-09-08, after a coordinator correction that was right. REV1 asserted "
            "PREFIX-hood ('the support is always an initial prefix of the energy order') and "
            "reported 0/2916 failures -- but every REV1 probability family carried an epsilon "
            "floor (+1e-30, 1e-18, +1e-6), so p>0 held everywhere and the assertion could not "
            "fire. REV2 adds exact-zero families, ASSERTS the true theorem (subset-hood + "
            "every hole is a zero-probability state) and MEASURES prefix-hood as a separate, "
            "regime-dependent statistic."),
        statement=(
            "The realised CVaR tail's support is always a SUBSET of an initial prefix of the "
            "energy order, and equals that prefix exactly when every state in the prefix "
            "carries positive probability. p_theta can delete a member; it can never add one "
            "outside the classical top-m."),
        n_cells=len(rows), n_failures=len(fails), failures=fails[:20],
        # the theorem
        subset_hood_violations=int(sum(1 for r in rows if not r["ok_subset"])),
        holes_that_were_not_zero_prob=int(sum(1 for r in rows
                                              if not r["ok_holes_are_zeros"])),
        # the over-claim, now a measured statistic rather than an assertion
        prefix_hood_violations=int(n_prefix_violations),
        prefix_hood_violation_rate=float(n_prefix_violations / max(1, len(rows))),
        # the empirical regime
        n_full_support_cells=int(n_full_support),
        n_full_support_cells_with_exact_equality=int(n_full_support_equal),
        distinct_tail_sets_per_realised_m=distinct_by_m,
        m_values_with_more_than_one_set=multi,
        n_random_p_draws=4000, random_p_draws_regime="full support (p > 0 everywhere)",
        verdict=("ALGORITHM: subset-hood is an identity of cvar_from_probs / tail_indices, "
                 "not a property of the pool; equality is the full-support regime"
                 ) if ok else "BROKEN",
        source_lines=dict(
            tail_indices="core/quantum.py:229-272 (k=ceil(alpha*n) line 252; mask from e alone)",
            cvar_from_probs=("core/quantum.py:336-359 (order 352, cum 353, take 354 -- the "
                             "clip's SECOND conjunct p[order]>0 is what punches holes)")),
        complete=True,
    )
    tmp = os.path.join(RESULTS, "d_setequality.json.tmp")
    with open(tmp, "w") as f:
        json.dump(out, f, indent=1, default=str)
    os.replace(tmp, os.path.join(RESULTS, "d_setequality.json"))
    print(json.dumps({k: v for k, v in out.items() if k != "failures"}, indent=1))
    return 0 if (not fails and not multi) else 1


if __name__ == "__main__":
    sys.exit(main())
