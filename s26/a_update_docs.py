"""s26/a_update_docs.py -- lane A: append the slide-qualifier table to agentA_FINDINGS.md,
refresh DELIVERABLES_CHECK.md, and add the STATUS lines. Idempotent on the findings section
(checks a sentinel)."""
from __future__ import annotations
import io, os, re, datetime
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
SENT = "Claims that reached a slide"

SEC = """

## Claims that reached a slide, and the qualifier the presenter says with each (as of {now}; slides per `s26/PRESENTATION_CHANGES.md`, L61 / L73 / L76)

| slide | claim on the slide | say with it | source of the qualifier |
|---|---|---|---|
| 2 | phi MAE 36.1 (full sequence context) vs 36.4 deg (sequence-blind) | pooled over 1,507 residues of 126 targets; the 0.3 deg gap is the whole measurable sequence-to-phi channel at this length; psi 62.4 vs 72.8 | `s26/results/a_c26_phi_mae.json`, L31 |
| 2 (notes) | the dev identity-leak price +0.0004 A | on the lam = 0 chain (`fit`), 0.33x its MDE; +0.0018 on the built chain with a CI spanning zero; the +0.0030 benchmark half is S10-4's historical figure and is not on the slide | L44, L50 |
| 3, 4 | 3.2148 A built chain (production) | the cache's mean; the results-lab rebuild reads 3.2126 (PDB round trip); quote one and name it | EXAMINATION H, L57 |
| 3, 4 | 3.0483 A raw average | a point cloud, 22.3% contracted, not a structure; never compared with a built chain | E section H, L28 |
| 4 | random-75 3.425, constant helix 4.065 | point-cloud and ladder bases; the built-chain twins of the physics rows are 3.2187 to 4.1015 | L28 item 3, L29 |
| 4 | benchmark +0.0103, CI [-0.160, +0.180] | the artefact is `s9/final_report.json`, never opened in S26; the two means beside it are pinned by passing tests; no validated improvement | C04-C06, L32 item 2 |
| 4 | the 2/60 self-copy bound 0.028 A | the mean-CI envelope under assumption A2; the same envelope's worst single target gives 0.151 A; MINOR under every reading; dev-proxy price 0.002 A; cannot move the benchmark verdict either way | L55, L58 |
| 5 | 0.902 nats from the Gibbs optimum; 5.6e-17; cos 1.000000000 | on the 78 alpha = 1 targets the optimum is a product state (KL to the product of marginals at most 7.9e-4, per-target rows); the readout is insensitive at 0.24x MDE | L68, L70 |
| 6, 9 | the S25 slopes (-0.649 to -0.047 log2 Var per qubit) and "no exponential plateau" | always "at depth 3"; the DLA is the full so(128) from depth 2, so the slopes are a shallowness statement, not a favourable-algebra statement; never "barren plateau" either way | R2, L45 |
| 6 | grown circuits give no width-scaling argument (A4) | at alpha = 1 the grown circuits are product circuits; the one non-trivial family (alpha = 0.25, L2) decays at -0.302 against the fixed -0.243 with no persisted slope CI, so "within error" is asserted | L47 |
| 7 | Legacy +0.330 / AMBER +0.455 worse than a random subset | point-cloud basis, both sides; 5/5 folds | L28, L29 |
| 7, 10 | "refine with physics is a validity step, not an accuracy step"; +0.0207 worse than doing nothing | +0.0207 at 2.16x MDE is clean; "+0.0111 worse than a random move" is Type-M (1.12x): say the sign, not the size; "validity step on 124 of 126" (2BP4 and 9KAR break a virtual bond; 9KAR does not converge) | L46 |
| 7 | AMBER as a steric reject filter is harmful (+0.228 at 1e4) | Type-M magnitude (1.17x); the direction is measured at 1e3 and by the shrink arm; the harm is tail-carried (median +0.003, worst +4.15); point cloud so far, the built chain pending | L54 |
| 7 (notes) | strain predicts error, Spearman +0.433 | a calibration flag, never a gain; "CIs replicated", not "result replicated" (a deterministic statistic on fixed rows); provisional until the pool-spread control lands | L81 |
| 8 | A1: ADAPT is null, -0.014 / -0.022 A at 0.23x / 0.36x MDE | the comparison resolves 0.06 A; "twelve of twelve arms negative" is one observation (pairwise correlation 0.955); the 21-parameter match is a budget with realised counts 7 to 21; growth at alpha = 1 is inert, not absent (L75) | L70, L75, L82 |
| 8 | the fixed circuit stops 0.90 nats short; ADAPT reaches the Gibbs state | on the 78 alpha = 1 targets; on the 48 alpha = 0.25 targets the Gibbs state is not the CVaR optimum (gibbs_T is +0.067 there) | L70 caveat 3 |
| 10 | the C2 anchor 3.2126 and the rung contrasts | rebuild basis (L57); noesm +0.208 is Type-M with the sign at 5/5; every other rung so far is under its MDE ("not measured", MDE beside it); any best rung at the close is an order statistic (`best_of_k_within`) | L79 |
| 11 | -- | no number of mine; the tie-break floor 0.024 A, if quoted, applies across runs that do not share the tie-break, not to paired within-run contrasts | L71 |

## OPEN (queue as of {now})

- The pool-spread control for L53 (`s26/a_strain_vs_spread.py`, job queued behind the four-job cap); L81 is provisional until it lands.
- Every new entry within the hour: the reject chain, A3, the remaining rungs (esm8m after L78's kill, mix, pairnet, raw), branch_select, rotamer_relief B, the ensembling, W's survivors, lane I's slow test tier and the results-lab rebuild.
- RETRACTIONS.md (R1-R5) and DELIVERABLES_CHECK.md kept current; the final Part 10 pass at about 04:00 with the rule-13 grep over the report, the deck dump and every presenter-facing file.
"""


def main():
    f = os.path.join(ROOT, "s26", "agentA_FINDINGS.md")
    t = io.open(f, encoding="utf-8").read()
    if SENT not in t:
        io.open(f, "a", encoding="utf-8", newline="\n").write(SEC.replace("{now}", now))
        print("findings: section appended")
    p = os.path.join(ROOT, "s26", "DELIVERABLES_CHECK.md")
    t = io.open(p, encoding="utf-8").read()
    t = re.sub(r"Last pass: [^\n]*", "Last pass: %s (after L82; sprint extended to ~04:30, L77)." % now, t, count=1)
    t = t.replace("`s26/RETRACTIONS.md` R1-R4 seeded, R2 marked FOR docs/FINDINGS.md",
                  "`s26/RETRACTIONS.md` R1-R5 (R5 = L75's inert-growth correction), R2 marked FOR docs/FINDINGS.md")
    t = t.replace("| deck / PRESENTATION_CHANGES | built (L61, slide 8 pending); grep at the final pass | -- | -- |",
                  "| deck / PRESENTATION_CHANGES | built (L61), slide 8 rebuilt (L73), L75 wording applied (L76); grep at the final pass | -- | -- |")
    io.open(p, "w", encoding="utf-8", newline="\n").write(t)
    q = os.path.join(ROOT, "s26", "STATUS.md")
    s = io.open(q, encoding="utf-8").read()
    head = "## A (Adversary)\n"
    i = s.index(head) + len(head)
    j = i
    while j < len(s) and s[j:j + 2] == "- ":
        j = s.index("\n", j) + 1
    line = ("- %s done: L79 (rungs L62/63/65/66/67/72 STAND; L66 selection wording caveat), L80 (L52), L81 (L53 provisional pending a_strain_vs_spread, queued at the job cap), L82 (my L70 caveat 2 corrected per L75); RETRACTIONS R5; slide-qualifier table in agentA_FINDINGS.\n"
            "- %s next: the spread control result into an L81 follow-up; then every new entry within the hour (reject chain, A3, esm8m/mix/pairnet/raw, branch_select, rotamer_relief B, ensembling, W survivors, I's slow tier and lab rebuild); final Part 10 pass ~04:00.\n") % (now, now)
    io.open(q, "w", encoding="utf-8", newline="\n").write(s[:j] + line + s[j:])
    print("deliverables and status updated")


if __name__ == "__main__":
    main()
