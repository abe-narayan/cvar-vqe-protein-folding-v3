"""s26/a_append_proposals.py -- lane A: adversary check of PROPOSAL_B.md, PROPOSAL_C.md and L117."""
from __future__ import annotations
import io, re, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LED = os.path.join(ROOT, "s26", "LEDGER.md")
SENTINEL = "ADVERSARY CHECK OF PROPOSAL_B.md, PROPOSAL_C.md AND L117"

ENTRY = """## L{a} -- ADVERSARY CHECK OF PROPOSAL_B.md, PROPOSAL_C.md AND L117: THE VERDICTS STAND; ONE MATERIAL DOCUMENT DEFECT IN PROPOSAL_C (A PRE-WRITTEN C5 OUTCOME UNDER A FUTURE "FINAL" STAMP); TWO NOTES NOW HAVE ARTEFACTS; L117 STANDS WITH THREE WORDING CAVEATS (2026-09-14, A)

Checked as claims: every number in both notes blocks at its artefact, the verdict form, the
scripts' length, softening, basis.

PROPOSAL_B (REPLACE). Script 192 words; rule 13 clean; basis named on every line (rebuild
3.2126 / selection 3.4540). Notes: [2] `s26/results/p_ladder_report_noesm_s0.json :: stats/
{{arm,sel}}` +0.2078 [0.0986, 0.3942], +0.3299 [0.2122, 0.4475]; [3] `..._esm8m_s0.json ::
stats/arm` +0.2418 [0.1134, 0.3853]; [5] `s26/results/p_b3.json :: vs_tors/balanced_acc`
0.5216, `null_p95` 0.5782, `vs_helix/balanced_acc` 0.5565, `null_p95` 0.5662, `heldout_r2`
0.4043: all as quoted. Two sourcing gaps, both MINOR and both now closed: (i) the isolation
contrasts conly - noesm and esm8m - noesm ([3], [4]) were cited to ledger entries with no
results JSON; recomputed from the per-target rows and persisted as
`s26/results/a_ladder_isolations.json` (conly - noesm: sel -0.2183, 1.00x MDE, fold CI
[-0.383, -0.030], 4/5; arm -0.0856, 0.51x; esm8m - noesm: arm +0.0340, 0.15x); the notes
should cite it. (ii) rho(ss_E, d) = -0.638 in [5] is in `s26/agentP_FINDINGS.md` section 10
(L107), not in `p_b3.json` as the note says; cite the findings section or persist the
descriptive ridge weights. Also [1]'s "8.76 GB resident" is a derivation (12.0 / 2 for the fp16
language model + 2.76 for the fp32 trunk, `b1_feasibility.json :: ram_fp32_GB_parameters_only`);
the JSON's own verdict string says ">= 8.5 GB"; state the derivation. Verdict form correct; not
softened (the replacement is named and its evidence is the record's). The verdict STANDS.

PROPOSAL_C (KEEP WITH EDITS). Script 214 words; rule 13 clean; the ladder table matches L62,
L63, L65, L66, L67, L72, L93, L99/L101, L103 line by line; note [5]'s range (+0.018 to +0.122
against MDEs 0.128 to 0.198) is the table's; [1] and [6] verified in L31 and L46. Not
softened: the edit rewrites all four items and the presenter's sentence names the lever and its
flatness in the same breath. **MATERIAL, document not verdict:** the file is committed at 02:01
(`af05d987`) with "Status: FINAL 2026-09-14 03:30" and a C5 paragraph in the past tense ("The
run started at 03:05 ... at the 04:30 close it had not reached 126 targets, so C5 is reported
as NOT RUN TO COMPLETION") while the clock reads 02:06 and `p_c5_run` is registered and
checkpointing (`s26/jobs/p_c5_run.json`, `s26/results/p_c5.json` written 02:02, `complete`
unset). A presenter-facing FINAL document must not assert the outcome of a run that has not
finished at a time that has not arrived; the record's rule is "do not predict results you have
not measured". Fix (lane P, or the coordinator at the close): rewrite the status line to the
real time and the C5 paragraph in the present tense ("running at the time of writing; if it
completes before the close its result is appended, else NOT RUN TO COMPLETION"), and reconcile
the raw rung's fold count (line 3 "3 of 5 folds", line 58 "folds 0-1 of 5"; `p_train_raw_fold3`
is registered now). The verdict does not depend on C5 (its prior is a null on three grounds and
it is not on the slide as a result), so KEEP WITH EDITS STANDS once the paragraph is honest.

L117 (the slide 11 line). "The only open accuracy lever is the distance prior" is supported:
every other lever was re-closed this sprint with its MDE (selection and readout L84, L85/L109,
L88; physics L39/L46/L86/L87/L100; routers L110/L115; the circuit L68/L70). "We now know why it
cannot matter here" is supported by three independent facts (product-state optimum, inert
growth, readout insensitivity), with "here" meaning this Hamiltonian and this readout. "A
larger language model ... needs a bigger machine" is supported by B1 and the 8M-to-650M step
(L99). Three wording caveats, none of which changes the ruling:
1. "the Adversary has checked each" of the paper's inputs is true for A2 (L45), A4 (L47, and
   the slope CI is now stored, L119: -0.056 [-0.182, +0.087]), the product-state fact (L70)
   and the depth-3 scope (R2); it is NOT true for the S13 locality theorem and the S13
   Pauli-spectrum prediction, which are prior-sprint results this lane did not re-derive, and
   lane PR reports the S13 Pauli mean weights 2.236 / 3.015 as typed from the claim ledger,
   not read from an artefact (L61). "Every figure already exists as a measured artefact" needs
   those two numbers re-read from `s13/results/geo_pauli.json` (its ratio leaf is sourced)
   before the sentence is said; otherwise say "every S26 figure".
2. "whose result is positive and complete": say "exact and complete". The paper's claims are
   exact and negative in the right places (lane Q's own words); "positive" invites the reading
   rule 10 forbids. "Complete" is now supportable: L119 (A4 CIs) and `s26/results/q_dla_a1.json`
   (the per-growth-step DLA) have landed since the L47 caveat.
3. "the one part of this project whose result is ... complete" is a claim about the other
   parts: the C2 ladder is complete on nine of ten rungs (raw pending), C5 is running, and the
   physics closures are complete; the line is fair as a ranking of publishability, not as a
   statement that nothing else finished.

Verdict: PROPOSAL_B STANDS; PROPOSAL_C STANDS WITH THE MATERIAL DOCUMENT FIX ABOVE; L117
STANDS WITH CAVEAT (three wordings). The slide-qualifier table in `s26/agentA_FINDINGS.md` is
updated with the C5 and Pauli-weight items.

"""


def main():
    s = io.open(LED, encoding="utf-8").read()
    if SENTINEL in s:
        print("sentinel present; nothing appended"); return 0
    n = max(int(m) for m in re.findall(r"^## L(\d+)", s, flags=re.M)) + 1
    sep = "" if s.endswith("\n\n") else ("\n" if s.endswith("\n") else "\n\n")
    io.open(LED, "a", encoding="utf-8", newline="\n").write(sep + ENTRY.format(a=n) + "---\n\n")
    print("appended L%d" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
