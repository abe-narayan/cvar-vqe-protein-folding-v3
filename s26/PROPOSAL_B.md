# PROPOSAL B -- A LEARNED FOLDING MODEL AS THE PRIOR (lane P, Sprint 26)

Status: FINAL 2026-09-14 01:15. B1, B2 and B3 measured. Verdict form per the campaign prompt.

## What the proposal says

Replace or augment the shipped distance prior with the output of a large pretrained folding
model (ESMFold), on the argument that the prior's accuracy is the only steep lever the project
has measured (S24 L13: -2.15 A per unit gamma at the origin) and that a model trained on the
whole PDB should be a better prior than a 250k-parameter MLP trained on 787 peptides. Its
repo-native fallback (B3) asks whether the targets on which the pipeline beats a sequence-only
predictor can be recognised without the native, so that the pipeline and a sequence-only model
could be routed.

## What the repository already knows

The prior's derivative is steep only along the native's own direction (S25 L12). The ESM-2 650M
representation is worth -0.288 A over one-hot on selection (S7-11; artefact lost, ledger L11,
re-measured this sprint at -0.330, L62) and its contact head is a filter, not a discriminator
(S17 L23). Five router constructions failed to route the per-target set size (S22 L7, S23 L7);
the S22 L10 bound puts a five-feature linear router at n = 281 for a 0.24 A effect. The pool's
error is 68% common-mode (S23 L9).

## B1 -- can the folding model run here? (final; `s26/results/b1_feasibility.json`, ledger L13)

No, on three independent grounds. (a) No ESMFold or ESM-2 3B checkpoint is on disk; the
download is 2,771,653,574 + 5,678,116,398 bytes (HTTP HEAD 2026-09-13; nothing downloaded).
(b) `import esm.esmfold.v1.pretrained` fails on `omegaconf` and `openfold` is absent; fair-esm's
README says openfold needs nvcc and Python 3.9 or lower; this box runs Python 3.13 with CPU-only
torch. (c) fair-esm loads the 3B language model in fp16 and the 690M trunk in fp32
(esmfold.py:43-46): 8.76 GB resident before the load transient against 4.4 GB of campaign
headroom (governor: 16.75 GB, 64.8% used, ceiling 93%); fp16 everywhere, 7.4 GB, does not fit
either. The 650M contact head is not a stand-in folding model: it is a contact map and is
already 13 of the shipped prior's 183 input columns.

## B2 -- the feasible-scale ladder (final; ledger L62, L63, L99/L101, `s26/results/p_ladder_report_{noesm,conly,esm8m}_s0.json`)

Every arm through one code path, paired against the shipped posterior; built chain on the
rebuild basis 3.2126 (L57), selection 3.4540. Negative = better than shipped.

| rung | inputs | built chain: effect (SE), MDE, fold CI, folds, W/L | selection | verdict (arm) |
|---|---|---|---|---|
| noesm | 42-d physicochemical pair block, no ESM | +0.208 (0.073), 0.205 (1.01x), [+0.099, +0.394], 5/5, 47W/79L | +0.330 (1.27x), [+0.212, +0.448], 5/5 | WORSE (Type-M zone) |
| conly | + the 13 contact-head columns, no embedding | +0.122 (0.071), 0.198 (0.62x), [-0.035, +0.267], 4/5, 53W/73L | +0.112 (0.42x) | NOT MEASURED |
| esm8m | esm2_t6_8M reps + its own contact head, for the 650M | +0.242 (0.074), 0.207 (1.17x), [+0.113, +0.385], 5/5, 52W/74L | +0.225 (0.88x), [+0.122, +0.305] | WORSE (Type-M zone) |
| pca32 | the 650M (shipped) | identity on 126/126 | identity | IDENTITY |

Isolations: conly - noesm (the contact head alone) -0.218 on selection (1.00x MDE, fold
[-0.374, -0.037], 4/5) and -0.086 on the built chain (0.51x); esm8m - noesm +0.034 (0.15x): the
8M model carries none of the channel. The size axis is not flat between 8M and 650M (-0.24 A on
the built chain), and the ladder cannot say whether 3B carries more; it can say that the value
is in the 650M representation, not in "any language model".

## B3 -- is the set where the pipeline gains over sequence-only characterisable native-free? (final; ledger L106/L107, `s26/results/p_b3.json`)

Built chains: pipeline 3.2148, sequence-only torsion predictor 3.7705, constant helix 4.0648;
tors - arm +0.556 (SE 0.132), 37W/89L; helix - arm +0.850 (SE 0.149). A nested leave-fold-out
ridge classifier of the SIGN of the gain on 45 native-free features is at its permutation null
against both comparators (balanced accuracy 0.522 vs a null 95th percentile of 0.578; 0.557 vs
0.566). The SIZE of the gain over the constant helix is partly predictable (held-out R2 0.40;
squared-error reduction 1.19x MDE, fold CI [-1.495, -0.700], 5/5) and the feature that carries
it is the retrieval pool's strand content (rho -0.64): when the pool is strand-like a constant
helix loses a lot, when it is helical it loses little. The sign is not predictable, so nothing
can be routed. FAIL18 (ORACLE stratum) is where the sequence-only predictor wins (+0.463) and
nothing native-free locates it. The pre-registered falsifier fires.

## Verdict: REPLACE

Per the campaign rule (coordinator, 2026-09-13): B1 cannot run on this box, B2 finds no
feasible-scale model input that beats the shipped prior, and B3 yields no native-free
characterisable set with a fold-clustered effect above MDE. The replacement is the
trainability paper outline written by lane Q, `s26/PROPOSAL_B_REPLACEMENT.md`. The
feasible-scale ladder facts (B2) fold into Proposal C as its second candidate: the ESM-2 650M
channel is the whole of the prior's learned input value this instrument can see (-0.208 A
built chain, -0.330 A selection, 5/5 folds), two thirds of its selection value is in the contact
head, and a larger language model is the one untested input, which needs a bigger machine.

## The two-minute script (presenter, sourced in the notes below)

"Proposal B wanted to swap our distance prior for a large pretrained folding model. We measured
whether that can even run here: it cannot. ESMFold needs about nine gigabytes resident and a
GPU-era dependency chain, and this machine has four gigabytes of headroom and no CUDA [1]. So we
measured the part of that idea that does fit. Our prior already reads a 650-million-parameter
protein language model; removing it costs 0.2 Angstroms on the built chain and 0.33 on
selection, five folds out of five [2]. Swapping in an 8-million-parameter model loses all of
that [3], and two thirds of the selection value sits in the model's contact head [4]. We also
asked whether we could tell, without the answer, which targets the pipeline beats a sequence-
only predictor on, so we could route between them. We cannot: a classifier is at chance, and
the only thing predictable is how much a constant helix loses, which is set by how strand-like
the retrieved pool is [5]. Our verdict is replace: the replacement is the trainability paper,
and the one input we could not test, a larger language model, needs a larger machine."

## Notes block

[1] `s26/results/b1_feasibility.json`; ledger L13: 8.76 GB resident vs 4.4 GB headroom;
`ModuleNotFoundError: omegaconf`, `openfold` absent; checkpoints 8.45 GB not on disk.
[2] Ledger L62, `s26/results/p_ladder_noesm_s0.json`: +0.208 (fold CI [+0.099, +0.394]),
+0.330 (fold CI [+0.212, +0.448]).
[3] Ledger L99/L101, `s26/results/p_ladder_esm8m_s0.json`: +0.242 (fold CI [+0.113, +0.385]);
esm8m - noesm +0.034 (0.15x MDE).
[4] Ledger L63: conly - noesm -0.218 on selection (fold CI [-0.374, -0.037]).
[5] Ledger L106/L107, `s26/results/p_b3.json`: balanced accuracy 0.522 / 0.557 vs null 95th pct
0.578 / 0.566; R2 0.404; rho(ss_E, d) = -0.638.
