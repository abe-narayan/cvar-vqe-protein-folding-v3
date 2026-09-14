"""s26/a_append_l84_l93.py -- lane A: adversary checks of L84-L91 and L93, four entries, atomic."""
from __future__ import annotations
import io, re, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LED = os.path.join(ROOT, "s26", "LEDGER.md")
SENTINEL = "ADVERSARY CHECK OF L84 AND L85"

ENTRIES = [
"""## L{a} -- ADVERSARY CHECK OF L84 AND L85 (W: window ensembling, window provenance): L84 STANDS; L85 STANDS WITH CAVEAT (AN ORACLE POOL FACT; THE DEPLOYABLE TEST H_P3 IS STILL TO RUN) (2026-09-14, A)

L84 (fixed-K ensembling of three BLOSUM keys). Pre-registered before the probe; basis stated
(built chain on the rebuild basis, point cloud carried, `sel` undefined for an ensemble and not
quoted); the zero-information control is matched in member count (three 75-resamples of the
production set); the primary is -0.0005 A at 0.02x MDE and the control +0.0042 at 0.14x, both
NOT MEASURED; every secondary is under its MDE, and the one that touches the Type-M edge (the
BLOSUM80 shortlist alone, +0.0237 at 0.70x with the fold CI above zero, 4/5) is flagged "if
anything worse", which is the right wording. Power: the design resolves 0.028 A and sees
nothing; the mandatory direction is closed at that resolution in its fixed-K form, the
widening-K form having been closed by S17 L12. Deterministic, nothing to replicate. STANDS.

L85 (window provenance). The census is native-free (codes and sequences); the contrast is
ORACLE on the single-window basis and labelled so. Pool: peptide-derived minus fragment -0.416
A, 3.18x MDE, fold CI [-0.447, -0.376], 101W/25L, 5/5, not concentrated: a clean ORACLE fact,
and the one S19 already implied (the peptide corpus carries the sequence-structure channel; the
fragments are 73% of every pool). Top-75: whole-or-terminal minus fragment -0.131 at 1.13x
(Type-M, flagged). Two caveats. (1) Nothing deployable is claimed yet: the class is native-free
but the contrast is ORACLE; the pre-registered ACHIEVABLE test H_P3 (class-weighted readout
against uniform and a permuted-weight control) is the only thing that could turn this into a
lever, and by the record (uniform readout optimal over two families, S25; 68% common-mode) the
expectation is null; the Adversary checks H_P3 when it lands. (2) Six class contrasts are
listed; the two the prereg named are the ones that count, and the four secondaries are under
their MDEs. STANDS WITH CAVEAT.

""",
"""## L{a} -- ADVERSARY CHECK OF L86 AND L87 (PH: the reject on the built chain; the C3 replication): BOTH STAND WITH CAVEAT; THE TOWARD-MEMBER GAIN IS A MEASURED, REPLICATED 0.02 A THAT UNDOES PART OF THE PROJECTION'S OWN MOVE (2026-09-14, A)

L86 (steric reject, built chain). The cross-basis replication L54 asked for: R@1e4 +0.248 A
against the re-projected shipped top-75 (1.17x MDE, Type-M, fold CI [+0.120, +0.348]) and
+0.157 against the same count rejected at random (1.15x, Type-M), S@1e4 +0.104 (1.11x, 5/5),
S vs RANDS +0.055 (0.73x, NOT MEASURED). Same caveats as L54, plus one: on the built chain the
R-arm primary holds on 4 of 5 folds (fold 1 -0.002), not 5/5, so it does not meet the standing
rule's fold clause on this basis; the S arm does (5/5), the 1e3 arms are measured on both bases
and the dose is monotone, so "closed on both bases, with a sign" stands on those. The harm is
again tail-carried (median +0.0016, worst +4.20 on 8T61). The anchor is the re-projected
shipped set (addendum 2), the right comparator for a re-projected arm. STANDS WITH CAVEAT.

L87 (C3 stage-1 replication). New seeds, reversed order, fold labels fixed: every replicated
quantity lands inside L39's fold CI. The toward-member control minus do-nothing is now -0.0207
[-0.0250, -0.0166], 5/5, 1.77x MDE, 94W/32L, replicated: by the discipline this IS a measured
result, and it is native-free (a move of AMBER's own magnitude toward a random member of the
shipped pool). PH reads it as a property of the projection's displacement, not of physics, and
not a proposal; I agree, with the mechanism stated in numbers so the presenter cannot mistake
it: the built chain sits 0.166 A further from the native than the point cloud it was projected
from (L38 / EXAMINATION H), the pool members surround that cloud, so a 0.22 A step toward a
member is a partial reversal of the projection's own displacement and recovers about an eighth
of it. The natural comparator, the same step toward the cloud itself, was not run and would
say whether "member" matters at all; and a moved chain is no longer an ideal-geometry chain,
so the emitted object changes basis. Not a lever; a measurement of what the projection costs.
The AMBER-minus-random Type-M status (1.02x in the replication) is unchanged. STANDS WITH CAVEAT.

""",
"""## L{a} -- ADVERSARY CHECK OF L88 (branch_select): STANDS; THE 56% SPLIT-HALF TRANSFER IS AGAINST THE COLUMN MEAN AND LANDS ON THE PRODUCTION CHOICE, NOT ABOVE IT (2026-09-14, A)

Operator native-free (the converged e1 of each of the five projection solutions; ties by
`ST.argmin_tied`, 0 ties reported); ORACLE evaluation only. The falsifier did not fire: e1 pick
minus production +0.0055 at 0.18x MDE with 42 exact ties (the e1 pick IS the production choice
on a third of the targets), the fold CI on the harmful side, 4/5 folds: NOT MEASURED, closed as
an accuracy step with the power stated. The e1 pick beats a random branch by -0.102 (2.29x,
5/5) and the raw single point picks worse (+0.092, Type-M): the relaxed energy carries the
same discriminating power among branches as the 2D torsion prior, and the raw energy carries
the builder's side-chain clash (L23). Clean on leakage, ties, both CIs, concentration and power.

One clarification, so the 56% is not read as a hidden lever. `ST.best_of_k_within` over the five
columns gives an ORACLE per-target minimum of -0.190 (3.130 A) against the row mean, with a
split-half transfer of -0.107 (56%, k_eff 4.63, argmin counts [33, 24, 39, 18, 12]). A
transfer is measured against the COLUMN MEAN, which here is the random pick (3.320 A); a fixed
best column therefore lands at about 3.21 A, which is the production objective's own choice
(3.213). So "a per-target-consistent column exists" means one start is systematically better
than a random start, and the objective already captures that; a fixed-start rule would at best
tie production. The ORACLE 0.083 A between production and the per-target best branch is the
order statistic L88 says it is. STANDS.

""",
"""## L{a} -- ADVERSARY CHECK OF L89, L90, L91, L93 (the cis floor's tight form; the recall gradient; memorisation on the ladder; the mix rung): ALL STAND; L38's 0.347 A IS NOW TO BE QUOTED AS THE OWN-TORSION UPPER BOUND BESIDE THE TIGHT 0.083 A (2026-09-14, A)

L89 (ORACLE DIAGNOSTIC). The native projected through the production projection sits 0.083 A
from itself (0.043 with the prior off; max 0.44); floor2 minus L38's own-torsion floor -0.264
(3.46x, 113W/13L); the prior costs +0.040 when the input is the native (2.39x, 23W/103L); the
chain cost exceeds floor2 by +0.083 (1.51x, 5/5). L38 declared its 0.347 an upper bound and
registered floor2 in advance, so nothing is retracted; every quotation of "0.35 A" now carries
"own-torsion upper bound; the tight floor is 0.08 A, and the projection's 0.166 A cost is the
operator's displacement, not representation" (added to the slide-qualifier table). STANDS.

L90 (native-free covariates against the ORACLE label). No covariate reaches the registered
|rho| >= 0.25; the strongest, max pinned identity to the corpus, is -0.205 (0.81x, fold CI
[-0.315, -0.099], permutation p 0.010; -0.176 after partialling n): the correlation analogue of
the Type-M zone, and W says so. The unnamed confound W adds (the training corpus IS the
retrieval library, so "recall" and "a nearer pool window" are not separable here) is the right
one and is stated before any reading. Power stated (|rho| >= 0.25 excluded at ~0.8). STANDS.

L91 (ORACLE DIAGNOSTIC, own-native models). The direction-discounted ladder (-2.1496 x gam x
cos) over-predicts the measured cloud delta by +0.346 A (1.50x MDE, 5/5, 65W/61L); the
registered falsifier fires, so the S25 L12 currency is calibrated on small re-readings only and
must not be used to price a large off-axis move. Not deployable by construction; the value is
the calibration, which L62 and L79 already apply ("gam_eff is never a prediction of the
endpoint"). STANDS.

L93 (mix rung). Leave-fold-out chooses lam = 0 on 5/5 folds, so the arm is the shipped
posterior on every basis (126 ties, the degenerate FLAG); every lam > 0 is monotonically worse
(+0.733 at lam = 1, typicality). The per-target ORACLE over the eight lams (-0.399) is 1.83x
covered by the valid across-target null and its 41% split-half transfer is the trivial statement
that lam = 0 beats the mean over lams; no lam other than 0 transfers. Nested by construction
(the lam chosen leave-fold-out), so no regression-to-the-mean concern. S7-3 / S19 L14's
foreshadowed null, measured with the falsifier. STANDS.

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
