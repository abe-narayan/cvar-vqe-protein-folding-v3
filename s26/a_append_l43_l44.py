"""s26/a_append_l43_l44.py -- lane A: append the adversary checks of L43 and L44 atomically."""
from __future__ import annotations
import io, re, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LED = os.path.join(ROOT, "s26", "LEDGER.md")
SENTINEL = "ADVERSARY CHECK OF L43"

ENTRIES = [
"""## L{a} -- ADVERSARY CHECK OF L43 (STERIC REJECT, POINT CLOUD): STANDS WITH CAVEAT (2026-09-13, A)

A harmful result, so the checklist is applied to the controls and the power, not to a gain.

- Leakage: the operators (`reject_refill`, `reject_shrink`, `random_shrink`, `random_refill`,
  `permuted_energy`, `retained_sets`, `s26/ph_reject.py:109-200`) read `e_amber`, `score_dist`,
  `sub` and `order` only; no `nat_ca` and no RMSD inside them (the `rr` at `:169` is an rng handle,
  not the ORACLE array; `univ_nativefree` refuses `rr`). The native enters only in `cloud_rmsd`,
  the ORACLE scoring of a native-free operator. Clean.
- Tie-breaking: the score order is `argsort_stable` (the production rule); an empty retained set
  falls back to the anchor and is counted as a tie (10 ties at 1e4 R, 13 at S), as pre-registered
  in addendum 1. Clean.
- iid vs fold CI: both quoted verbatim; at 1e4 both exclude zero for R vs anchor and R vs RANDR,
  5/5 folds. Clean.
- Concentration / median-vs-mean: the uniform-effect null is computed and not flagged (pctile
  0.52). But the harm is TAIL-CARRIED: R@1e4 vs anchor has median +0.0027 against mean +0.2276,
  49W/67L/10T, p90 +1.28, worst +4.15 (8T61). The presentation must say "near zero on the median
  target, catastrophic on a minority (the 8 targets whose whole pool has no survivor, and the deep
  refills)", never "+0.228 on every target". CAVEAT 1.
- k_eff for the threshold sweep: `ST.best_of_k_within` applied (k_eff 3.74 R / 3.72 S), and the
  negative split-half transfer is correctly read as "the mildest threshold transfers", not a gain.
  Clean.
- Tuned parameter: the primary 1e4 was fixed before the run (PREREG section 3). Clean.
- Baseline / basis: the anchor is the production `rmsd_avg` (reproduced to 1e-6 on every
  target); controls matched in count, refill and permutation; point cloud on both sides, stated.
  Clean.
- Type-M: the two HEADLINE numbers, R@1e4 vs anchor +0.228 (1.17x MDE) and R vs RANDR +0.167
  (1.26x), are Type-M-zone magnitudes. The direction "harmful" rests on the measured contrasts:
  1e3 R +0.532 (1.75x), 1e3 S +0.196 (1.66x), 1e4 S +0.108 (1.34x), RANDS +0.037 (1.47x), 5/5
  folds throughout, and on the monotone dose. Quote the 1e4 R magnitudes with the flag. CAVEAT 2.
- "Refilling from ranks 76 to 147 is nearly free" rests on RANDR vs anchor +0.061 at 0.57x MDE
  (UNDERPOWERED, iid CI includes zero), so the refill-cost / choice-cost decomposition is a
  point-estimate split, not a measured one. CAVEAT 3.
- Replication: a negative; the cross-basis replication is the built-chain job `ph_reject_chain`,
  running; the direction must hold there. Power stated correctly for the 1e5 / 1e6 nulls.

Verdict: STANDS WITH CAVEAT (the conclusion -- the physical-threshold reject with refill is
harmful on the point cloud, its limit is the anchor -- is established; the +0.228 / +0.167
magnitudes are Type-M, the harm is tail-carried, and the refill-cost split is underpowered).
Disposition "closed on the point cloud in the two unmeasured forms" is accepted pending the
built chain.

""",
"""## L{a} -- ADVERSARY CHECK OF L44 (THE 2/60 PROXY BOUND): STANDS WITH CAVEAT; THE CLASS MINOR IS STABLE UNDER EVERY READING OF THE ENVELOPE, THE NUMBER 0.028 IS NOT (2026-09-13, A)

- Leakage / reads: the native-free halves use `load_blind` (`rr`, `nat_ca` NaN-poisoned); the
  gated `endpoint` and `floor` read natives for scoring only; no benchmark file, sequence or name
  anywhere (confirmed in L49 and re-checked on `s26/w_endpoint_report.py`'s inputs). Clean.
- Ties: `emit` uses `argsort(kind="stable")`, the argmin tie set is stored and `sel` is averaged
  over it (`ST.argmin_tied`'s rule); the refill follows the production corpus order. Clean.
- iid vs fold: Part C on arm -0.698, fold CI [-0.831, -0.576], 2.83x MDE, 112W/14L, 5/5;
  concentration at the null's 52nd percentile; the median (-0.266) sits well inside the mean
  (-0.698), a broad but skewed effect. Clean, ORACLE-labelled, and L50 rightly keeps it out of
  the price of the leak that exists.
- Bound arithmetic: every row of `s26/results/w_selfcopy_bound.json` recomputes: (2/60) x
  0.0683 = 0.0023 (2P5H); (2/60) x 0.2016 = 0.0067 (triangle); (2/60) x 0.8312 = 0.0277 (arm
  envelope fold-CI limit); (2/60) x 1.4377 = 0.0479 (sel); (2/60) x 0.6838 = 0.0228 (gain).
  Materiality thresholds pre-registered (0.017 / 0.170). Clean.

Two caveats.

1. **Artefact / ledger disagreement on the paired gain.** `w_selfcopy_bound.json :: verdict/gain`
   reads IMMATERIAL at 0.0006 (the n = 1 `both_removed4` row) because `C_envelope_fold_ci` carries
   no `gain` row and `report()`'s `max(cands)` (`w_selfcopy.py:913`) therefore never saw the
   envelope for that basis; L44's Part D table and its pre-registered verdict say MINOR at 0.023
   (from the gain row +0.5101 [+0.3698, +0.6838]). The ledger's class is the prereg-correct one
   (Part D: IMMATERIAL only if B_real AND B_env are both below 0.017). Lane W should add the
   envelope gain row to the artefact so the report cites a JSON that agrees with the ledger.
2. **B_env is the fold-CI limit of a MEAN effect, not a per-target bound.** The prereg defined it
   so (Part D, assumption A2), and A2 is stated, but a bound on the contribution of two SPECIFIC
   targets is (2/60) x a per-target quantity. From `s26/results/w_selfcopy_endpoint.json :: C/rows`
   (leaked minus clean, arm, mean over the four leaked models): median -0.266, p05 -2.490, worst
   -4.536 A (2BP4); over (target, model) pairs the worst is -4.557. So the same envelope gives
   (2/60) x 4.536 = **0.151 A** as the worst-single-target bound and (2/60) x p95 = 0.083 A;
   the artefact's own native-free triangle rows say the same (`C_envelope_p95_over_targets/arm`
   0.147, `C_envelope_max_over_models/arm` 0.211, the latter loose because branch flips move
   the chain orthogonally to the native, as L44 notes). The class MINOR therefore holds under
   every reading of the envelope (0.028 mean-CI, 0.083 p95, 0.151 worst target; all below
   0.170), which is a stronger statement than L44 makes; but "bounded at 0.028 A" must be quoted
   as "expected contribution 0.023, mean-CI limit 0.028, worst single target 0.151, under A2",
   and the envelope named as the over-bound it is (own-native training, 60x the one measured
   carrier effect).

Also noted, not a caveat: Part B is 1 of 10 retrains (host kill, L40), stated as a deviation;
F2's "at least half" missed by one (5 of 11) and is reported as such; F1, F3 (n = 1), F4, F5
hold. Verdict: STANDS WITH CAVEAT. The benchmark caveat text should read: "2/60 self-copies;
dev-proxy price 0.002 A; own-native envelope 0.028 A (mean CI) to 0.151 A (worst target), MINOR
under every reading; cannot move the benchmark verdict either way."

""",
]


def main():
    s = io.open(LED, encoding="utf-8").read()
    if SENTINEL in s:
        print("sentinel present; nothing appended"); return 0
    n = max(int(m) for m in re.findall(r"^## L(\d+)", s, flags=re.M))
    block = ""
    for e in ENTRIES:
        n += 1
        block += e.format(a=n) + "---\n\n"
    sep = "" if s.endswith("\n\n") else ("\n" if s.endswith("\n") else "\n\n")
    io.open(LED, "a", encoding="utf-8", newline="\n").write(sep + block)
    print("appended L%d through L%d" % (n - len(ENTRIES) + 1, n))
    return 0


if __name__ == "__main__":
    sys.exit(main())
