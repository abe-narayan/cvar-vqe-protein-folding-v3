"""s26/a_append_rungs.py -- lane A: append four adversary entries atomically (rungs; L52; L53; L70 correction)."""
from __future__ import annotations
import io, re, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LED = os.path.join(ROOT, "s26", "LEDGER.md")
SENTINEL = "ADVERSARY CHECK OF THE C2 RUNGS"

ENTRIES = [
"""## L{a} -- ADVERSARY CHECK OF THE C2 RUNGS L62, L63, L65, L66, L67, L72: ALL STAND; L66's "A REAL, SMALL DEGRADATION" ON SELECTION IS NOT LICENSED AT 0.48x MDE; THE LADDER'S EVENTUAL BEST RUNG IS AN ORDER STATISTIC (2026-09-13, A)

Common path checked once (`s26/p_ladder.py`): the phase gate is enforced in code (`:118`
`SIGNOFF`; `:64`); native quantities enter only in the read-after-selection half (`:513-522`,
"ORACLE-READ-AFTER-SELECTION"), after `score_target`; the top-75 is `argsort(kind="stable")`
(`:505`); `sel` is tie-averaged with `ST.argmin_tied` (`:517`); every rung goes through one code
path against the one anchor (`p_ladder_shipped_s0.json`), whose built chain is the rebuild basis
3.2126 (L57) and is bit-exact against the cache on selection and point cloud; gate 2 (L65) shows
the retrained pca32 emits the shipped answer on 126/126, so rung differences are inputs or
architecture, not trainer noise. iid and fold CIs on every contrast; concentration nulls
computed; medians beside means; power stated. Clean on the checklist.

Per rung:
- L62 noesm: WORSE at 5/5 folds on all three bases, but every magnitude is Type-M (1.01x /
  1.02x / 1.27x): the sign is measured, the "+0.208 A" and "+0.330 A" are upper bounds; L62
  says so. The harm is tail-carried (median +0.044, worst +3.65 on 8T61); the uniform-effect null
  does not flag it. "Reproducing S7-11's -0.288 to -0.34" is a reproduction of a Type-M number
  by a Type-M number; quote both with the flag. The gam_eff-positive-while-worse reading (S25
  L12 from the training side) is correct and important: gam_eff is never a prediction of the
  endpoint.
- L63 conly: correctly UNDERPOWERED (0.62x), not null; the "two thirds of the ESM channel is in
  the contact head" is a point-estimate split (-0.218 at exactly 1.00x MDE against -0.330 at
  1.27x, both Type-M); the built-chain split is unresolved and L63 says so. The worst target is
  1CEK (+3.70), the verbatim self-copy target: noted, not read into.
- L65 pca32: 126 exact ties on all bases; the degenerate all-zero concentration FLAG is the
  identity case. The strongest possible gate-2 pass.
- L66 wide: built chain UNDERPOWERED (0.25x). CAVEAT: the selection sentence "a real, small
  degradation, worth reporting as 'the wider head selects slightly worse'" is not licensed:
  +0.097 is 0.48x its MDE, the iid CI spans zero ([-0.051, +0.240]) and only the fold CI excludes
  it; by the standing rule (clear the MDE AND fold CI excluding zero) it is not a result and the
  word is "suggestive, not measured". The conclusion "more capacity does not buy a better prior"
  stands on the built-chain null and the flat MAE.
- L67 pca32f: every contrast at 0.11x to 0.17x; "changes nothing measurable" is correct at the
  stated MDEs (0.13 to 0.18 A) and closes the global-PCA leak question at that resolution.
- L72 pca128: 0.43x / 0.38x / 0.13x, all under MDE; the isolate against pca32f (+0.058, 0.37x)
  is the clean size-axis point; power stated (a 0.18 A gain would be seen).

Multiplicity, stated now for the ladder's close: seven rungs so far and four to come, all
paired against one anchor on one pool, so the per-rung deltas are correlated (the rungs share
the pool, the score and the projection); any "best rung" or per-target minimum over rungs must
go through `ST.best_of_k_within` with its k_eff and split-half transfer (PREREG_C2's own rule),
never the raw minimum. No rung is positive so far, so nothing is owed yet. Verdicts: L62, L63,
L65, L67, L72 STAND; L66 STANDS WITH CAVEAT (the selection-basis wording).

""",
"""## L{a} -- ADVERSARY CHECK OF L52 (conformational_identity_floor, ORACLE DIAGNOSTIC): STANDS WITH CAVEAT (n = 18; BOTH PAIRED CONTRASTS ARE TYPE-M) (2026-09-13, A)

Reads: the copy's CA-RMSD to the target's model-1 native at the shared segment and the pool's
ORACLE `rr`; every quantity reads the native and is labelled ORACLE DIAGNOSTIC; nothing selects;
dev targets only (lane I's L15 list, 18 targets, 22 partners; no benchmark sequence, confirmed
in L49 for the same script). Registered falsifier F5 (median above 1.5 A; FALSIFIED below 1.0)
holds at 2.908 A over 22 pairs (18% below 1.0 A, 27% below 1.5 A); the four cross-fold
self-copies reproduce S24 L4's four values. Deterministic (Kabsch on fixed coordinates; the
bootstrap is seeded). Power: n = 18 gives MDEs of 0.92 and 1.26 A on the two paired contrasts.

Caveat. Both paired contrasts are in the Type-M zone (copy minus pool best +0.968 at 1.06x MDE;
copy minus pool mean -1.443 at 1.14x); W labels them so ("sign clean, magnitude not a result").
The number for the report is the median 2.9 A, a descriptive statistic at n = 22, not the paired
effects; and the sentence "the 2.0 A target is below what a verbatim lookup reaches" is a
statement about this instrument's 18 targets with a verbatim relative, not a general law. W's
own "not done" (no length-matched non-verbatim control population) is the right missing arm;
the pool mean plays that role loosely. Verdict: STANDS WITH CAVEAT.

""",
"""## L{a} -- ADVERSARY CHECK OF L53 (strain_difficulty): STANDS WITH CAVEAT, PROVISIONAL UNTIL THE POOL-SPREAD CONTROL LANDS; "REPLICATED" MEANS THE CIs, NOT THE ESTIMATE (2026-09-13, A)

Checked against `s26/ph_strain.py` and `s26/results/ph_strain.json` / `ph_strain_rep.json`.
- Leakage: the four signals are the production cache's own relaxation scalars (`amber_moved`,
  `amber_e0`, `amber_e1`, `amber_strain_after`), native-free by construction; `rmsd_arm` enters
  only as the ORACLE label (`:84-90`); gated (`:74`); nothing is selected, weighted or tuned.
  Clean.
- Multiplicity: four pre-registered signals, Bonferroni (p < 0.0125); the falsifier fires on one
  (`moved`, partial rho +0.433, fold CI [+0.247, +0.588], 5/5 folds, permutation p < 0.00025).
  Confounds n and Rg partialled by residualisation; per-fold signs and a 4000-draw permutation
  null. Clean. The two 4/5-fold signals (log_e0, log_drop at +0.24, +0.25) are correctly reported
  as measured-but-failing-the-rule, and `strain_after` as null; the FAIL18 Fisher tests as null.
- Outliers: rho +0.41 without 9KAR and 2BP4 (the broken-bond emissions). Clean.
- Basis: built chain, ORACLE label, stated. Clean.

Two caveats.
1. **"REPLICATED" in the title means the bootstrap and permutation draws, not the estimate.**
   The Spearman is a deterministic function of 126 fixed rows; the registered replication
   (new seed, reversed order) can only re-draw the CIs and the permutation null, and it did
   (identical rho +0.433, fold CI [+0.248, +0.581]). It cannot test whether the signal
   generalises to new targets; that needs targets the record does not have (no fresh
   benchmark exists). Say "CIs replicated", not "result replicated".
2. **The prereg carried no pool-disagreement control.** The obvious native-free difficulty
   proxy is the top-75's own pairwise CA-RMSD spread (the members' disagreement, no native).
   If `moved` is that spread in disguise, L53's "first native-free quantity above 0.4" is a
   re-labelling. The Adversary's control `s26/a_strain_vs_spread.py` (rho of the spread with
   `rmsd_arm`; partial rho of `moved` given n, Rg AND the spread; job `a_strain_vs_spread`,
   queued behind the four-job cap at 22:3x) decides it; the verdict below is provisional until
   it lands and is confirmed or amended in a follow-up entry.

Verdict: STANDS WITH CAVEAT, provisional. As a calibration flag with the L53 wording ("a
calibration curve, never a gain") it may be quoted now; the novelty sentence waits for the
control.

""",
"""## L{a} -- CORRECTION TO THE ADVERSARY'S L70 CAVEAT 2 PER L75: L-BFGS DOES APPEND OPERATORS AT alpha = 1; THEY ARE INERT; THE CAVEAT'S CONCLUSION IS UNCHANGED (2026-09-13, A)

L70 caveat 2 said "on the 78 alpha = 1 targets L-BFGS grows nothing and the 'P21' arm is the
7-parameter RY layer". Per L75 (lane Q, from the records' `adapt.*_lbfgs_zrank.sequence`), under
L-BFGS operators are appended on 60 of 78 (pool V) and 68 of 78 (pool L2) targets, all
multi-qubit, worth at most 1.2e-4 nats, with angles at most 0.018 rad, leaving the state a
product to 4.1e-4 nats. My `len(ops)` reading (0 to 21 under L-BFGS) was consistent with that
and I mis-stated it as "grows nothing". Corrected wording of caveat 2: the parameter match is a
budget match (21) and the realised counts run 7 to 21 under both optimisers, with the L-BFGS
appended operators inert; nothing in the A1 null depends on the count (P7, P14, P21 agree within
0.01 A). Recorded in `s26/RETRACTIONS.md` R5 together with lane Q's own corrections.

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
