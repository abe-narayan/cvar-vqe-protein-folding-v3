# Sprint 31 — CVaR-VQE Protein Folding: the objective, the readout, and where the information is not

**Status: DRAFT IN PROGRESS.** Sections marked `[PENDING]` await lanes still running. Assembled as
results land so nothing is reconstructed from memory at the end.

Branch `s26` · instrument: `tuning126`, 126 targets, 9–16 aa · endpoint: **mean built-chain Cα RMSD**
Charter: `s31/BRIEF.md` (verbatim, 1,660 lines) · Ledger: `s31/LEDGER.md` · State: `s31/STATE.md`
Contract: `s31/S31_CONTRACT.md` (32 rules) · Verifier: `s31/s31_verify.py` · Multiplicity:
`s31/MULTIPLICITY.md`

---

## 0. The answer, up front

[PENDING — written last.]

---

## 1. Endpoint and reproducibility

### 1.1 The endpoint did not move

**Production remains 3.2105 Å** (`s29/results/s29_O_chain_rows.jsonl :: item=prod`), mean built-chain
Cα RMSD over the 126 `tuning126` targets. The CA point cloud is **3.0483 Å** and the *set mean* is
**3.5507 Å** — three different objects, and this report names the basis in the same sentence as every
number.

`benchmark60` is **sealed**; its single pre-registered pass is spent and two guard sites enforce it
(`core/pipeline.py:386-388`, `core/bench.py:1074`, `core/pipeline.py:1664`). Folds and clusters were
not regenerated.

### 1.2 The projection is deterministic and chaotic — there is no seed, and the record said otherwise for three sprints

S28, S30's ledger, S30's report, S30's verifier and this sprint's own charter all describe *"an
unpinned multi-start projection seed."* **There is no RNG anywhere on the projection path.**
`core.project.fit_multi` loops over four fixed starts and takes a strict argmin; `fit_prior` is a
deterministic L-BFGS-B call. Reprojecting the same cloud twice is **bit-identical — max |ΔCA|
exactly 0.0**, unchanged under BLAS thread count.

**The real defect is conditioning, and it is worse than a seed would be:**

```
1A13's four lam=0 objectives   0.5091924865 / 0.5091925170 / 0.5091924488 / 0.5091924188
                               -> spread 1e-7, i.e. NOISE
the same four branches at lam=0.3   0.947 / 0.519 / 0.562 / 2.396
                               -> spread 1e-1
```

A decision taken where the objective **cannot** discriminate fixes an outcome where it **can**.
Measured across the benchmark: the median λ=0 branch margin is **3.1e-7**, and **73 of 126 targets
(58%)** have their branch chosen below a 1e-6 margin. A **1e-14 relative** perturbation of the input
cloud moves the emitted chain by a median of **1.6e-3 Å with a 0.511 Å tail** (2LNG), reproducing the
historical 0.517 Å discrepancy exactly. **Amplification ~1e13. Not one target of 126 is unchanged to
1e-9.**

**Three rules follow and were binding on every lane from the moment they were measured:**

1. Reprojection is reproducible **only from bit-identical clouds**. Both sides of any built-chain
   contrast must be projected **in the same job from the same stored clouds**.
2. **Canonical stays 3.2105 Å** — the run the endpoint was declared from, and the only one whose
   clouds are persisted per target. Not retro-fitted.
3. **No built-chain claim below 0.0107 Å** — the spread across the five circulating values.

**A consequence nobody had drawn:** the measured 0.92 cloud→chain transfer is **an average over a map
that is locally chaotic on a substantial minority of targets.** A lane proposing a cloud-level gain
must expect a non-smooth chain response there.

### 1.3 The reference itself is uncertain — and it does not explain the tail

**92.9% of `tuning126` is NMR-determined** (115 solution NMR, 5 electron crystallography, 4 X-ray, 2
solid-state NMR; RCSB GraphQL over all 126). The manifest reference is *"deposited coordinates, model
1"* — an arbitrary member of an ensemble. The ensembles are genuinely wide: mean pairwise CA-RMSD
between deposited models **1.0823 Å**, with model 1 sitting **0.6965 Å** from its own ensemble medoid
on average (max 4.294); 55/111 resolved targets have spread > 1.0 Å and 15/111 > 2.0 Å.

**It does not explain the tail:**

```
corr(production RMSD, ensemble spread)  +0.1118   95% CI [-0.0762,+0.2921]   SPANS ZERO
worst 18    mean RMSD 6.0291   mean ensemble spread 0.9672
other 93    mean RMSD 2.6842   mean ensemble spread 1.1046
tail-minus-rest  -0.1374   SE 0.2421   MDE 0.6783  ->  0.20x MDE   NULL
```

**The tail's targets have, if anything, narrower deposited ensembles.** The same null holds in all
four other arms. FAIL18 is real failure against a reference no worse determined than any other
target's.

**What survives is a caveat on the absolute number only:** a **uniform ~0.70 Å reference term** — a
perfect predictor aiming at the ensemble medoid still scores ~0.70 Å against model 1. In quadrature
`sqrt(3.21² − 0.70²) = 3.13` against 3.21, i.e. **~0.08 Å today and material only near 1 Å**. Because
it is uniform it **cancels in every arm-to-arm delta**, which is this project's actual currency.
**This is explicitly not a reason to re-score against the medoid** — the reference is part of the
sealed instrument. Reported as an uncertainty line; nothing changed.

---

## 2. Scientific hypothesis

[PENDING]

---

## 3. Mathematical mechanism

[PENDING]

---

## 4. CVaR-VQE architecture

[PENDING]

---

## 5. Hamiltonian

[PENDING]

---

## 6. Quantum encoding

[PENDING]

---

## 7. Objective geometry

[PENDING]

---

## 8. Free-energy work

[PENDING]

---

## 9. Torsion work

[PENDING]

---

## 10. Sparse / readout work

[PENDING]

---

## 11. Candidate-index information

[PENDING]

---

## 12. FAIL18 analysis

[PENDING]

---

## 13. Controls and nulls

[PENDING]

---

## 14. The 126-target endpoint

[PENDING]

---

## 15. Statistical analysis

[PENDING]

---

## 16. What was closed

[PENDING]

---

## 17. What improved

[PENDING]

---

## 18. What did not improve

[PENDING]

---

## 19. ORACLE ceilings

[PENDING]

---

## 20. The remaining information bottleneck

[PENDING]

---

## 21. Next-sprint recommendation

[PENDING]

---

## Appendix A — every claim withdrawn this sprint

[PENDING]

## Appendix B — the five "worst 18" strata, as a key

Five different strata are called "the worst 18" somewhere in this project's record. They differ by up
to **0.81 Å on the same-named quantity**. Every such number in this report names its stratum in the
same sentence.

| stratum | defined by | `best1_500` |
|---|---|---|
| `FAIL18` | the filter's own recall (`s12/instrument.py:271-278`) | 2.2842 |
| `defn18` | top-18 by widening gain | **= FAIL18, 18 of 18 — an identity** |
| `worst18_poolmean` | pool mean | 2.2286 |
| `worst18_bestpool` | best pool member | 3.0930 |
| "the genuinely worst 18" | production built-chain RMSD | ORACLE best pool member 2.5298 |
| worst 18 by production RMSD, n = 111 matched | ensemble arm | mean RMSD 6.0291 |
