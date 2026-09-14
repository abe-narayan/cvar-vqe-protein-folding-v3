"""s26/a_append_l64_l68.py -- lane A: append the adversary checks of L68 and L64 atomically."""
from __future__ import annotations
import io, re, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LED = os.path.join(ROOT, "s26", "LEDGER.md")
SENTINEL = "ADVERSARY CHECK OF L68"

ENTRIES = [
"""## L{a} -- ADVERSARY CHECK OF L68 (A1, qubit-ADAPT vs THE FIXED ANSATZ): STANDS WITH CAVEAT; THE REPLACE VERDICT (L69) IS SUPPORTED (2026-09-13, A)

Checked against `s26/results/a1/<pdb>.json` (126 records), `s26/results/a1_stats.json`,
`s26/q_adapt.py` and `s26/PREREG_A1.md`.

Confirmed clean.
- Readout: every arm goes through the production functions `pl.consensus_medoid`,
  `pl.average_weighted` and `pl.project` on its own `p` (`s26/q_adapt.py:849-851`); the
  comparator `fixed_zrank_it50` reproduces the S25 quantum-arm cache bit-for-bit (`cache_check`:
  `ca` and `q_ca` max abs 0.0, `sel_equal` True, 126 of 126). Basis built chain both sides.
- Leakage: the native enters only in `label_one` (`:887-911`, RMSD scoring), gated on "PHASE 0
  SIGNED OFF" (`:925`); the build half reads no native. Clean.
- Registered threshold: "ADAPT is null fires if both PRIMARY contrasts lie within +-0.5x their
  own MDE" (`s26/PREREG_A1.md:85`); measured 0.23x and 0.36x. Fires as registered.
- Power statement: complete and correct (MDE 0.0588 / 0.0608, Gelman-Carlin power 0.10 / 0.18,
  Type-M 3.7 / 2.4, "null at the registered threshold; underpowered below 0.06 A"). Fold CIs
  span zero, 3 of 5 folds, concentration at the null's 50th percentile. No replication owed.
- The product-state claim rests on persisted per-target rows: `product_diagnostics/zrank/
  kl_gibbs_to_product` (recomputed over 126: mean 1.40e-4, max 7.90e-4, as quoted) and
  `arms/fixed_zrank_it50/kl_to_gibbs` on the 78 alpha = 1 targets (mean 0.9027, max 0.9836, as
  quoted; ADAPT-L2-P21 2.6e-4 mean, 7.0e-4 max). Sourced.

Three caveats.
1. **"Twelve of twelve arms negative" is one observation, not twelve.** The twelve ADAPT delta
   vectors (built chain, per target) have mean pairwise correlation 0.955 (min 0.919, max 1.000)
   and their first principal component carries 96.0% of their variance; `ST.best_of_k_within`
   over the 13 arms (fixed + 12) gives a per-target oracle of -0.094 A whose split-half transfer
   is +0.006 (-6% of the oracle: NOT A SIGNAL). L68 rightly calls none of them a result; the
   phrase "a consistent direction" must be read with the 0.955 beside it, so that twelve
   correlated draws at 0.2-0.4x MDE are not taken for a trend.
2. **The parameter match is a budget match, not a realised one.** `max_params = 3n = 21` for the
   arms labelled P21, but `len(ops)` runs 14 to 21 for Adam and 0 to 21 for L-BFGS (`adapt/*/ops`
   in the records): on the 78 alpha = 1 targets L-BFGS grows nothing and the "P21" arm is the
   7-parameter RY layer, and Adam grows 7 to 14 strings. Nothing in the null depends on it (P7,
   P14 and P21 agree within 0.01 A on both bases), but "21 parameters both sides" should be
   stated as "the same 21-parameter budget; realised counts 7 to 21".
3. **The `gibbs_T` control's +0.0088 overall masks two opposite subsets.** On the 78 alpha = 1
   targets the exact Gibbs state is -0.0271 against the fixed circuit and ADAPT-L2-P21 is
   -0.0223; the two agree with each other (+0.0047, SE 0.0060, 0.28x MDE, fold CI [-0.0009,
   +0.0098]) as two states 2.6e-4 nats apart should. On the 48 alpha = 0.25 targets `gibbs_T` is
   +0.0671 (the Gibbs state is not the CVaR optimum there) and ADAPT is +0.0000. So the
   product-state diagnosis is internally consistent on the subset it applies to; and, noted for
   the record, 29 of 78 alpha = 1 targets differ by more than 0.01 A between two states 2.6e-4
   nats apart (max 0.206 A): the projection readout amplifies sub-milli-nat differences on some
   targets, the mechanism L64 measures as the pipeline's own noise floor.

Verdict: STANDS WITH CAVEAT. The REPLACE verdict for Proposal A (L69) is supported: the
endpoint is null at the registered threshold with a stated resolution of 0.06 A, the optimum is
a product state on every real target by persisted per-target rows, and no arm family changes
the emitted structure beyond a third of the MDE.

""",
"""## L{a} -- ADVERSARY CHECK OF L64 (THE TIE-BREAK NOISE FLOOR): STANDS WITH CAVEAT (2026-09-13, A)

- Operator and reads: the boundary tie class is re-drawn uniformly with `s15.seed.stable_rng`
  (8 seeds), the non-tied prefix untouched, everything downstream production; the draws are
  native-free (`load_blind`) and the endpoint is gated; production gate top-75 == `sub` on
  126 of 126. Clean.
- Numbers: m_tie 0.0039 A (sd of the 126-mean over 8 draws; W states its ~27% relative SE),
  paired MDE between two draws median 0.0236 A over the 28 pairs (max 0.0321), per-target s_tie
  median 0.0227, p90 0.129, above 0.1 A on 21 of 126; 3.57 of 75 members replaced per draw; the
  argmin unchanged on 91.5% of cells. Registered predictions held (m_tie, paired MDE) or were
  over-estimates (8 of 75 replaced; median s_tie 0.03-0.08). `s26/results/w_tiebreak_report.json`.
- The production convention against a random draw: 0.15x MDE on the built chain, 0.55x on the
  cloud, 0.75x on `sel` with 107 ties and the fold CI excluding zero -- correctly NOT MEASURED,
  correctly read as suggestive of S25's non-neutral retrieval order and not a lever.
- Replication: the 8 seeds are the replication; no fit, no fold order. Power: the floor is the
  power statement.

Caveat. The title's sentence "EVERY HUNDREDTHS-LEVEL EFFECT ON THE RECORD IS INSIDE THE LATTER"
is true only as W's body qualifies it: the floor is the noise between two RUNS that do not share
the tie-break (a corpus re-order, a re-retrieval, a cross-instrument or cross-sprint comparison).
Every one of the nine listed effects (C3's +0.0207 and +0.0111, S24's -0.022, the 0.015 ORACLE
weight, the +0.004 reranking, C27's two prices, L19's 0.012, L24's +0.013) was measured as a
PAIRED contrast with the pool and its tie-break held fixed, and their paired SEs (C3: 0.0034)
already exclude this noise; none of them is invalidated. The presenter must quote the report
sentence W wrote ("a hundredths-level effect is real only as a paired contrast with the
tie-break held fixed; quoted across runs or against another instrument it is inside the
pipeline's own convention noise, 0.024 A at n = 126"), never the title alone. Second, minor:
m_tie from 8 draws carries a 27% relative SE, so "0.004" is 0.003 to 0.005; the verdict does
not depend on it.

Verdict: STANDS WITH CAVEAT (the title needs the body's qualifier wherever it is quoted).

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
