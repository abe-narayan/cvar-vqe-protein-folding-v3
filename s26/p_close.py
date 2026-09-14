"""One-off close-out: C5 ledger entry (ST.fmt verbatim), raw not-run entry, PROPOSAL_C addendum 2, findings closing sections."""
import json, datetime, re, os, sys
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, '.')
from s24 import stats_lib as ST
now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
st = json.load(open('s26/results/p_c5_stats.json'))
def blk(a, b):
    o = dict(st['%s|%s' % (a, b)]); o['per_fold'] = {int(k): v for k, v in o['per_fold'].items()}
    return ST.fmt(o)
arms = ['global_R1_a0.5', 'global_R2_a0.5', 'ridge_R1', 'random_R1', 'random_R2_a0.5', 'oracle_R1', 'oracle_R2']
text = "\n\n".join(blk(a, 'arm') for a in arms) + "\n\n" + "\n\n".join(blk(a, 'cloud') for a in ('global_R1_a0.5', 'global_R2_a0.5', 'ridge_R1', 'random_R2_a0.5', 'oracle_R2'))
L = 's26/LEDGER.md'; t = open(L, encoding='utf-8').read(); n = max(int(m) for m in re.findall(r'^## L(\d+)', t, re.M)) + 1
c5 = ("\n## L%d -- C5 COMPLETE: PREDICTING THE COMMON-MODE DIRECTION ON HELD-OUT FOLDS AND SUBTRACTING IT IS NULL (GLOBAL) TO HARMFUL (RIDGE +0.164 A, 1.79x MDE, 5/5 FOLDS); THE ORACLE CEILING IS -1.87 A (DISTANCE SPACE) AND -3.13 A (COORDINATE FRAME); THE RIDGE PREDICTOR IS WORSE THAN A RANDOM MOVE OF ITS OWN SIZE (%s, lane P)\n\n"
      "`s26/p_c5.py run`, `s26/results/p_c5.json` (complete: true, 126/126; job `p_c5_run`, exit 0, 5,284 s, peak RSS 0.09 GB),\n"
      "statistics `s26/results/p_c5_stats.json`; pre-registered in `s26/PREREG_C5.md` (alpha grid reduced to 0.5, declared in the\n"
      "commit of `p_c5.py`). Construction: per target the shipped top-75 cloud c; the common mode ebar = c - t after one rigid\n"
      "Kabsch fit (ORACLE, training folds only inside the fits). R1 = distance space, five SHELL means of D(c) - D(t), applied\n"
      "by stress descent initialised at c; R2 = the cloud's own principal-axis frame (native-free; round trip exact, rotation-\n"
      "invariant to 3e-15 in the selftest). GLOBAL = the training-fold mean correction per length band; RIDGE = nested\n"
      "leave-fold-out ridge from the 45 B3 features to the five R1 shell means (alpha chosen 100 to 10,000 per fold, i.e.\n"
      "near-constant predictions); RANDOM = a random direction of the same magnitude through the same operator; ORACLE =\n"
      "the true correction (the bound). Every arm through `s12.instrument.project`; built chain primary (rebuild basis\n"
      "3.2126), point cloud carried. Negative = better than the incumbent.\n\n```\n%s\n```\n\n"
      "Reading. (1) The ceiling is large: subtracting the TRUE common mode gives 1.34 A (distance space, -1.87, 120W/6L)\n"
      "and 0.08 A (the coordinate frame, -3.13, 126W/0L; the R2 oracle is the native by construction, so its number is\n"
      "the size of the mode plus the projection's cost, not an achievable target). (2) The two achievable arms are null:\n"
      "GLOBAL R2 -0.009 (0.23x MDE, 3/5 folds, 70W/56L) and GLOBAL R1 +0.031 (0.60x, fold CI above zero at 5/5 but under\n"
      "the MDE gate). (3) The learned arm is HARMFUL: RIDGE R1 +0.164 A [fold +0.118, +0.225], 1.79x MDE, 5/5 folds,\n"
      "40W/86L, and its magnitude-matched random control is +0.135 [+0.029, +0.246], 1.33x, so the predicted shell profile\n"
      "is indistinguishable from a random shell profile of the same size (+0.029, inside the noise). The nested alphas\n"
      "went to the top of the grid on 3 of 5 folds: the ridge found nothing to fit and emitted a near-constant\n"
      "correction, which is the GLOBAL arm plus noise. (4) Power: the design resolves 0.04 to 0.05 A on the GLOBAL arms\n"
      "and 0.09 A on the ridge arm; a real gain of 0.05 A would have shown as a fold CI below zero on GLOBAL R2 and did not.\n"
      "Verdict: the pre-registered falsifier fires for both GLOBAL and RIDGE; C5 is refuted at one agent-day, as S16,\n"
      "S19 L14 (\"the estimator is made of the bias\") and S24 L7 predicted, now with the ceiling measured beside it on the\n"
      "same operator: the common mode is 1.9 to 3.1 A of built-chain error, and nothing native-free in these two\n"
      "representations touches it.\n\n---\n" % (n, now, text))
raw = ("\n## L%d -- THE raw RUNG (1280-d ESM-2 WITH A LEARNED PROJECTION) IS NOT RUN: 4 OF 5 FOLDS TRAINED, FOLD 4 UNTRAINED, NO EVALUATION; MEMORY AND TIME (%s, lane P)\n\n"
       "`s26/models/p_ladder/raw_fold{0,1,2,3}_s0.pt` exist (jobs `p_train_raw_f1`, `p_train_raw_fold1/2/3`: 4,477 to 4,900 s\n"
       "per fold at peak RSS 1.93 to 1.96 GB; the first attempt was terminated at 19:11 with the governor dead, L42, and\n"
       "fold 3 was queued 40 min behind the job cap). At the 04:30 close fold 4 is untrained and the rung cannot be\n"
       "evaluated (a rung needs all five leave-fold-out models). No further launch: a fifth fold is 80 min of a 2 GB job\n"
       "and the evaluation another 10, past the close. What the record says instead: raw's inputs are bracketed by pca32f\n"
       "(32 components, +0.018 built chain, 0.13x MDE, L67) and pca128 (128 components, +0.076, 0.43x, L72), both null,\n"
       "and S7-11 measured raw on selection as indistinguishable from pca32 and pca128. The rung is recorded as NOT RUN,\n"
       "not as a result; its four checkpoints stay on disk for the next sprint.\n\n---\n" % (n + 1, now))
open(L, 'a', encoding='utf-8').write(c5 + raw)
# PROPOSAL_C addendum 2
p = 's26/PROPOSAL_C.md'; s = open(p, encoding='utf-8').read()
s = s.replace("Status: last edited 2026-09-14 02:12 (corrected per ledger L120; see addendum 1).", "Status: last edited %s (addendum 2: C5 final)." % now)
s += ("\n## ADDENDUM 2 (%s) -- C5 FINAL: REFUTED; RAW NOT RUN\n\n"
      "`s26/results/p_c5.json` (complete: true, 126/126) and `s26/results/p_c5_stats.json`; ledger L%d. Subtracting a\n"
      "predicted common-mode correction on held-out folds: GLOBAL (coordinate frame) -0.009 A on the built chain, 0.23x MDE,\n"
      "3/5 folds, NOT MEASURED; GLOBAL (distance space) +0.031, 0.60x, NOT MEASURED; RIDGE (distance space, nested\n"
      "leave-fold-out on 45 native-free features) +0.164 A, fold CI [+0.118, +0.225], 1.79x MDE, 5/5 folds, WORSE, and\n"
      "indistinguishable from a random correction of the same size (+0.135, 1.33x). ORACLE ceilings -1.87 A (distance\n"
      "space) and -3.13 A (coordinate frame): the common mode is most of the error and nothing native-free touches it\n"
      "(S16, S19 L14, S24 L7 confirmed on the same operator). Edit 3 of the verdict therefore reads: item C5 is CLOSED at\n"
      "one agent-day, its ceiling measured. The raw rung (1280-d input) is NOT RUN (4 of 5 folds trained; ledger L%d);\n"
      "its inputs are bracketed by the null pca32f and pca128 rungs.\n" % (now, n, n + 1))
open(p, 'w', encoding='utf-8').write(s)
# findings closing sections
with open('s26/agentP_FINDINGS.md', 'a', encoding='utf-8') as fh:
    fh.write("\n## 12. C5 result and the close (ledger L%d, L%d; %s)\n\n"
             "DEMONSTRATED (D12): predicting the common-mode direction native-free and subtracting it is null (GLOBAL R2\n"
             "-0.009, 0.23x MDE; GLOBAL R1 +0.031, 0.60x) to harmful (RIDGE R1 +0.164, 1.79x MDE, 5/5 folds, 40W/86L; its\n"
             "magnitude-matched random control +0.135, 1.33x). ORACLE DIAGNOSTIC: the true common mode subtracted is\n"
             "worth -1.87 A (distance space, 120W/6L) and -3.13 A (coordinate frame, 126W/0L) on the built chain. Artefacts\n"
             "`s26/results/p_c5.json`, `s26/results/p_c5_stats.json`. raw: NOT RUN (4 of 5 folds trained).\n\n"
             "## 13. What damaged my own expectations (closing)\n\n"
             "- I expected the sign of gam_eff to track the sign of the endpoint. It does not: every null-to-worse rung has\n"
             "  positive gam_eff (+0.10 to +0.13, cos 0.18 to 0.29), and PairNet's cos 0.53 with an 11% MAE cut buys +0.041.\n"
             "- I expected the contact head to be worth little as a prior input (S17 L23 called it a filter); it carries\n"
             "  two thirds of the ESM channel's selection value (-0.218 vs noesm) and the 8M model carries none of it.\n"
             "- I expected the retrained recipe to differ from the pinned models by noise; it emits the pipeline's answer on\n"
             "  126/126 targets, which made the ladder cleaner than planned.\n"
             "- I expected the arm-choice oracle over {pipeline, sequence-only, helix} to carry per-target signal on FAIL18;\n"
             "  the across-target null accounts for 92% of it (L14), and B3's classifier is at chance.\n"
             "- I expected my own clock to be right. It ran 80 minutes fast for two hours and put a future stamp on a\n"
             "  proposal (L120/L124); the fix is in the record.\n\n"
             "## 14. What I did not do and why (closing)\n\n"
             "- raw rung: not evaluated (fold 4 untrained; 80-minute folds at 2 GB behind a 6-7 job cap; ledger L%d).\n"
             "- Coherence rungs (PREREG_coherence, p_coh.py, synthetic-tested): not run; no slot before the close.\n"
             "- Replication (second seed, reversed fold order): never triggered; no rung cleared its MDE in the improving\n"
             "  direction with 5/5 folds.\n"
             "- `attn` (IDEA_better_prior_inputs): the box never emptied; the 650M forward pass with head weights was not\n"
             "  probed.\n"
             "- C5 alpha grid reduced from three values to one (0.5) for the close, declared before the run.\n"
             "- No benchmark was opened; no pinned file was written; `esm_cache.npz` was loaded once, under jobrun, as a\n"
             "  probe (L11).\n" % (n, n + 1, now, n + 1))
print("posted L%d (C5) and L%d (raw not run); PROPOSAL_C addendum 2; findings 12-14" % (n, n + 1))
