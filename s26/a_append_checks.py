"""s26/a_append_checks.py -- lane A: append the adversary-check ledger entries atomically.

Re-reads s26/LEDGER.md, finds the current max L number, and appends the five check entries
(L27; L39; L35; L22-L24 + L38; L30) as consecutive `## L<n>` headings in one write, so no
other lane can interleave within the block. Run once; it is idempotent only in the sense that
re-running would append duplicates, so it checks for a sentinel first.
"""
from __future__ import annotations
import io, re, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LED = os.path.join(ROOT, "s26", "LEDGER.md")

SENTINEL = "ADVERSARY CHECK OF L27"

ENTRIES = [
"""## L{a} -- ADVERSARY CHECK OF L27 (A2, THE DLA): STANDS (2026-09-13, A)

Independent re-derivation `s26/a_dla_check.py` -> `s26/results/a_dla_check.json` (88 s, run
directly, under 200 MB; written from the construction in `s26/PREREG_A2.md` and
`core/quantum.py`, sharing no closure code with `s26/q_dla.py`, and using a different numeric
rank rule -- incremental Gram-Schmidt, no shared svd tolerance). Symbolic closure at n = 4..7,
L = 1..4 equals `s26/results/q_dla.json :: results/fixed/*` at every cell: n7 dim(DLA) = 8128 =
so(128) from L = 2 (7 at L = 1); n6 510 / 1023 / 2016; n4 120 from L = 2; n5 496 from L = 2. All
four combinations of conjugation convention (right, left) and gate order (forward, reverse) give
8128 at n = 7, L = 2 and L = 3. The independent numeric route agrees with the symbolic count at
n = 4, 5, L = 1..3 (6 of 6); every closure is odd-Y. The 8128 / 1025 reconciliation in
`s26/agentQ_FINDINGS.md` 2.1 is correct: 8128 is `q_dla.json :: results/fixed/n7_L3/dim` (the
fixed ansatz, so(128), 1.0 of so); 1025 is `results/adapt_sets/adapt_L2_*_n7_a0.25_T0.3/
ladder[P=21]/dim` (the strings ADAPT selected), 1025 / 8128 = 0.126. Leakage / ties / iid-vs-fold
CI: not applicable (an exact count, no native, no RMSD). Q's own H2b prediction (dim < 8128,
guess 4095) is correctly recorded as FALSIFIED. Verdict: STANDS. The "S25 slopes are a statement
about depth 3" reading is an inference from Larocca 2022 / Ragone 2024, not a measurement here;
it is a scope note for RETRACTIONS, not a retracted number.

""",
"""## L{a} -- ADVERSARY CHECK OF L39 (C3 STAGE 1): STANDS WITH CAVEAT (2026-09-13, A)

The three contrasts reproduce from `s26/results/ph_c3_stage1.json` (complete, 126 rows,
provenance e480fc15): AMBER minus do-nothing +0.02069 (effect/MDE 2.16, fold CI [0.0154,
0.0290], 5/5, WORSE), AMBER minus matched RANDOM +0.01108 (effect/MDE 1.12, fold CI [0.0062,
0.0171], 5/5), AMBER minus TOWARD-MEMBER +0.03851 (effect/MDE 2.27, fold CI [0.0298, 0.0479],
5/5); cos(AMBER, true residual) mean -0.0491 (SE 0.0149, 36.5% positive). Checks:

1. Control construction against S16. `s26/ph_lib.py:random_displacement` reproduces
   `s16/repair.py`'s `rand` exactly: an isotropic Gaussian on the CA trace, the six rigid-body
   directions removed by projection onto `s15/align_lib.py:rigid_basis`, scaled so
   ||g||/sqrt(n) equals AMBER's per-atom RMS displacement. The magnitude matched is `mag_sup`
   0.2198 A (SE 0.0084), measured after superposition, the same definition S16 used. Confirmed.
2. The Type-M reading. AMBER-minus-RANDOM is at 1.12x MDE with `type_m_flag` True (0.7-1.3x is
   the Type-M zone). So the HEADLINE +0.0111 is a Type-M number: its sign and fold CI are clean
   (5/5, CI excludes zero) but its magnitude is inflated ~1.07x and by the discipline it is
   "not a result" as a magnitude. The robust statements are AMBER-minus-do-nothing (+0.0207,
   2.16x) and AMBER-minus-member (+0.0385, 2.27x), both clear of the Type-M zone. Any
   presentation line must quote +0.0111 with the Type-M flag or lean on the two clean contrasts.
3. The decision rule holds regardless: "accuracy step" required AMBER to BEAT both controls;
   it is worse than both, so it is not an accuracy step. Note that TOWARD-MEMBER minus
   do-nothing is -0.0178 (BETTER, fold CI [-0.0255, -0.0124], 5/5): a zero-information move
   toward a random pool member IMPROVES RMSD where AMBER's physics move worsens it.
4. "Validity step" wording. Sound for the 124 targets that converge with a sane virtual bond,
   but on 2BP4 (relaxed CA-CA 5.38 A) and 9KAR (4.86 A, e1 1262 > CONVERGE_MAX_KCAL, not
   converged) the relaxation BREAKS a virtual bond; on those two it is not cleanly a validity
   step either. The presentation should say "a validity step on 124 of 126; on 2 it breaks a
   virtual bond, one of which does not converge."

Verdict: STANDS WITH CAVEAT (the +0.0111 headline is Type-M; "validity step" carries the 2BP4 /
9KAR exception). The finding's direction and its central conclusion are correct.

""",
"""## L{a} -- ADVERSARY CHECK OF L35 (A4, GROWN-CIRCUIT VARIANCE): STANDS WITH CAVEAT (2026-09-13, A)

The reproduction gate holds: `s26/results/q_var.json :: results/reproduction/passed` True,
`worst_rel` 0.0 (the S25 n = 7 rows of all five cells reproduce bit-identically). The slopes
reproduce (`results/slopes`): fixed deployed_a1_T03 -0.311, grown V +0.024, grown L2 -0.008;
the one non-trivial family, deployed_a025_T03, fixed -0.243, grown V -0.239, grown L2 -0.302.
Checks:

1. "Product circuit at alpha = 1" is established, but by `s26/results/q_dla.json`, not by
   `q_var.json`. `q_var.json` carries only the variance slopes; the product-circuit reading
   rests on the ADAPT closures in `q_dla.json :: results/adapt_sets` (alpha = 1 selects an
   abelian, dim-7 = n set, `n_distinct_ops` 1 to 2, the ladder dim stays 7), which I
   independently reproduced in `s26/results/a_dla_check.json` (L27 above). A large gradient
   from a product circuit is the trivial regime, not trainability (rule 10's mirror). Sound,
   with that provenance noted.
2. The one non-trivial slope's comparison. L35 says grown L2 -0.302 is "the fixed ansatz's
   -0.243 within the sampling error of a 7-point slope (relative SE of a variance 9 to 16%)."
   `q_var.json` stores the point slope but no CI on it, so "within error" is asserted, not
   computed: the difference is 0.059 in log2-slope-per-qubit and no SE on that difference is on
   disk. The qualitative conclusion is safe (both slopes are order -0.25 to -0.30, both far from
   the -1.0 of a 2-design), but the precise "within error" claim is not backed by a stored CI.
   Also: the ledger table's grown-V value -0.246 for that cell is the P = 3n-only 6-point fit,
   while `q_var.json :: slopes/deployed_a025_T03/grown_V_adam_best` is -0.239; the two differ
   because of the degenerate-row exclusion the ledger states, not a discrepancy.

Verdict: STANDS WITH CAVEAT (the negative conclusion -- ADAPT gives Proposal A no width-scaling
argument -- holds; the "within error" of -0.302 vs -0.243 is asserted without a persisted slope
CI, and "product circuit" is grounded in q_dla.json, not q_var.json).

""",
"""## L{a} -- ADVERSARY CHECK OF L22, L23, L24, L38 (PH CENSUSES AND THE CIS FLOOR): ALL STAND; L22 WITH A STALE-LINE-NUMBER CAVEAT (2026-09-13, A)

L22 (cis census, `s26/results/ph_cis_census.json`, complete). Omega statistics recomputed from
its rows: 1,507 bonds (= sum of n-1), mean 1.91, median 0.34, p90 5.8, p99 17.4, max 42.8 deg
at 9UV5 (bonds -137.2 and +144.8), 0.53% beyond 20 deg, 0.13% beyond 30; 0 of 126 model-1
natives cis by either criterion (the two agree 126/126), 0 of 1,966 ensemble models, 0 of
2,352,893 windows, universe minimum step 3.5045 A. Native-free path confirmed (natives through
`ph_lib.native_backbone`, omega/CA only; bank through `univ_nativefree`, which refuses `rr` /
`nat_ca`; record through `prod_record_nativefree`, RMSD keys stripped). CAVEAT: L22 cites the
step gate at `core/data.py:406-407` and `697-698`; at HEAD it is `core/data.py:439`
(`step.min() < 3.5 or step.max() > 4.1`) and `:725` -- the numbering before lane I's `37bddbbb`.
Same code, stale line numbers. Q's registered ensemble-cis prediction is correctly recorded as
FALSIFIED (0 of 1,966). Verdict: STANDS WITH CAVEAT.

L23 (steric reject census, `s26/results/ph_reject_census.json`, complete). The pool-identity
assertion is live code (`s26/ph_reject.py:96-102`: `universe_idx == I.pool_idx`,
`amber_verify_max_rel == 0.0`, production `sub` == score top-75 as a set) and holds on 1A13,
1S9Z, 9KAR (`s24/cache_amber/<pdb>.npz :: universe_idx` equals `order[:500]`; top-75 members
above 1e4 = 12 / 27 / 74, equal to the census rows; pool frac above 1e4 = 0.390 / 0.650 /
0.888). 96.8% = (2,627 + 2,267) / 5,057 side-chain-involving; rho(e, min heavy-atom distance)
within top-75 -0.743 (SE 0.017). Verdict: STANDS.

L24 (C3 native-free part). Every number recomputed from `bench_results/cache/1fc9f2dcf489e2fb`
with my own Kabsch: displacement 0.2198 A (SE 0.0084, median 0.1972, range 0.103-0.591),
`amber_moved` 0.2332, e0 median 8.59e4 / min -473 / max 1.3e14 / 58.73% above 1e4, e1 mean
-559.8 (SE 30.0) / max 1262.4 (9KAR), 125 of 126 converged, strain 60.4, relaxed bond mean
3.867 (min 3.12 at 1M02, max 5.38 at 2BP4, 4.86 at 9KAR), 6 of 126 outside [3.6, 4.0], Rg
+0.0457. All match L24. Verdict: STANDS.

L38 (cis floor, ORACLE DIAGNOSTIC, `s26/results/ph_cis_floor.json`, complete). floor_ca mean
0.3468 A (SE 0.029, median 0.272, max 1.474 at 1ID6); rebuild_bb 0.3400; chain_cost 0.1664;
rho(floor_ca, max omega dev) +0.828, rho(chain_cost, omega dev) -0.036, rho(chain_cost,
floor_ca) +0.083; the paired contrast cost-minus-floor -0.1804 (fold CI [-0.2266, -0.1022],
verdict BETTER). Every quantity reads the native and is ORACLE-labelled; the "BETTER" is
explicitly disarmed in L38 (it means only that the production chain cost is smaller than the
own-torsion rebuild floor). floor_ca is declared an UPPER bound on the manifold floor and the
tight `floor2` is registered. The two Spearman nulls (+0.08, -0.04) exclude |rho| above ~0.25
at n = 126. Verdict: STANDS.

""",
"""## L{a} -- ADVERSARY CHECK OF L30 (W's 2/60 PROXY-BOUND PREREG): SOUND, STANDS (2026-09-13, A)

L30 is a pre-registration plus a native-free census, not an endpoint result; the check is on
its soundness and on whether the census reads anything it must not.

Reads. `s26/w_selfcopy.py census` iterates `I.targets()` (the 126 dev targets), loads every
universe through `load_blind` (which overwrites `rr` and `nat_ca` with NaN, lines 105-115), and
reads the production record through `I.shipped_record` (native-free). No benchmark sequence,
name, PDB, native, RMSD or manifest is read anywhere in the census, retrieval, envelope or
posterior commands; the only native reads (`I.load_univ` at lines 686, 821) are inside the
gated `endpoint` and `floor`, which refuse to run before "PHASE 0 SIGNED OFF". The one
benchmark-derived fact used, 2/60, is on the record (S24 L4, lane I L15/L18). Confirmed clean.

Assumption set (A1-A5). Reasonably complete for the quantity claimed (the leak's contribution
to a benchmark mean, bounded in absolute value). A1 (mechanism match, no interaction when one
carrier carries both benchmark targets) and A5 (the envelope bounds channel B in the HELP
direction only; HARM is measured on the dev 4 directly) are the two load-bearing assumptions and
both are stated. The gap a reader would press -- that the benchmark carriers might be MORE
homologous to their targets than the dev-4 carriers, which would make the benchmark effect
larger than "max of four" -- is covered by the ORACLE-insertion arm (8 dev targets with a
same-fold verbatim carrier, longer identity 0.6 to 0.93 against the dev 4's <= 0.52), reported
beside the dev 4. The honest limitations are self-declared: n = 4 supports no quantile ("max of
four" is named as such), and the AMBER stage's second-order contribution to a per-target delta
is an assumption (A3), not a measurement.

Verdict: STANDS as a sound pre-registration. The bound's headline will rest on n = 4 for the
realistic arm and on the n = 126 envelope for the guard; when the gated endpoints land, the
Adversary re-checks the signed deltas, the triangle bounds and the envelope's fold CI against
this prereg before any benchmark caveat text is written.

""",
]


def main():
    s = io.open(LED, encoding="utf-8").read()
    if SENTINEL in s:
        print("sentinel present; entries already appended, doing nothing")
        return 0
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
