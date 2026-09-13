# IDEA (lane Q) -- PUBLISH THE TRAINABILITY HALF (open item 3 of the state brief)

## Hypothesis

The S13 and S25 trainability measurements, with A2 (the DLA) and A4 (grown-circuit variance)
from this sprint, form a paper whose claims are all exact, cross-validated and negative in the
right places: the energy model does not reshape the variational manifold; its Pauli spectrum
predicts gradient variance with no free parameter through a kernel that is flat in Pauli
weight at the deployed depth; the CVaR non-linearity flattens the width decay; the deployed
ansatz's Gibbs target is a product state; and the readout cannot tell a trained state from
its optimum 45% of its mass away.

## Why the record does not already close it

The record holds every measurement but no manuscript. `docs/CONDENSED_REPORT.md` ("What
survives: the trainability result") and `s25/QUANTUM.md` are internal. Open item 3 of
`docs/STATE_BRIEF_2026-09-12.md` names publication as open, and S25 QUANTUM.md section 7.1
names the DLA as the one unmeasured object a referee would ask for first. A2 measures it.

## Exact falsifier

This is a writing task, not an experiment; its falsifier is editorial. The outline
(`s26/TRAINABILITY_PAPER_OUTLINE.md`) is "publishable as a workshop contribution" only if every
figure has an artefact path, every claim a scope condition, and the two forbidden claims are
absent. The coordinator or examiner can fail it on any figure without an artefact.

## Expected effect

None on RMSD. The deliverable is `s26/TRAINABILITY_PAPER_OUTLINE.md` with a figure list, the
venue-honest gap list (hardware absent; noiseless exact simulation; n <= 13; no 2-design;
conditioning of raw force fields), the two claims not to make, and the related work.

## Memory and agent-hours

No compute beyond A2 and A4. Outline 2 agent-hours (done in this turn); a full draft 2 to 3
agent-days, outside this sprint.
