# Sprint 32 — Where the RMSD is lost, and what it would take to get it back

**Status: DRAFT IN PROGRESS.** Sections marked `[PENDING]` await lanes still running. Assembled as
results land so nothing is reconstructed from memory at the end.

Branch `s26` · instrument `tuning126`, 126 targets, 9–16 aa · endpoint **mean built-chain Cα RMSD**
Charter `s32/BRIEF.md` (verbatim, 2,087 lines) · Ledger `s32/LEDGER.md` · State `s32/STATE.md`
Contract `s32/S32_CONTRACT.md` · Causal map `s32/CAUSAL_MAP.md` · Verifier `s32/s32_verify.py`
Multiplicity `s32/MULTIPLICITY.md` · Theory `s32/THEORY_Q.md`

---

## 0. The answer, up front

[PENDING — written last.]

---

## 1. What the endpoint is, exactly

The charter's Step 4 asked for the 3.2105 Å endpoint to be verified from artefacts rather than taken
from a report. It was, and the answer is more precise than the number.

**3.2105 Å** is *the λ = 0.3 multi-start projection arm of `s12/instrument.project` applied to
`s29/results/s29_O_structs/<pdb>.npz["prod"]`, CA-RMSD to native, meaned over `tuning126`.*
**It is not a cached scalar.** Five distinct objects live near it and they are **not five estimates
of one**:

| value | what it is |
|---|---|
| 3.048338 | CA point cloud — the top-75 coordinate average, unprojected (`rmsd_avg`) |
| 3.204076 | the **λ = 0** arm — nearest ideal geometry, no Ramachandran penalty (`rmsd_fit`) |
| 3.214765 | the **λ = 0.3** arm — **the chain production itself emits** (`rmsd_arm` / `ca`) |
| 3.235460 | that chain **after AMBER relaxation** (`rmsd_full`) — AMBER costs **+0.0207 Å** |
| **3.210534** | **the canonical endpoint** — a *re-projection* of the stored cloud, λ = 0.3 |

**The canonical endpoint is 0.0043 Å better than the chain production emits.** A reader who assumes
the endpoint is the pipeline's output is wrong by more than several historical claims are large.

### 1.1 Everything upstream of the chain reproduces bit-exactly

Targets 126; folds 25/23/25/23/30; the top-75 **set** reproduced 126/126 from an independent
rescoring; the recomputed coordinate average matches the stored one to **1.42e-14**; cloud
**3.048338**; set mean **3.550683**; pool best K=500 **1.710824**; top-75 best **2.306153**; shipped
argmin **3.454000** to six decimals.

### 1.2 The projection operator is bit-reproducible; its *input* is not uniquely determined

| input cloud | mean chain | mean `d` | max \|d\| | bit-identical |
|---|---|---|---|---|
| `s29_O_structs/<pdb>.npz['prod']` | 3.210533995 | +0.000000 | 0.000000 | **126/126** |
| production cache `avg_ca` | 3.214765154 | +0.004231 | 0.416765 | 0/126 |
| recomputed coordinate average | 3.212625220 | +0.002091 | 0.517410 | 0/126 |

All three are **the same top-75 coordinate average, agreeing to 5.7e-14**, with cloud RMSDs identical
to 9 dp. `s12.instrument.project` reproduces **bit-for-bit across processes, across BLAS thread
counts 1/2/4/8, and across the gap between S29 and today**. Lane R reproduced five separate ladder
rungs bit-for-bit from a different job, script and process.

> **The operator is bit-reproducible. The input is what is fragile.** A **one-ULP (7.1e-15 Å)**
> change in the cloud moves the built chain by **0.10–0.15 Å** — deterministic given identical bits,
> **discontinuous** in them.

**Consequence:** the anchor carries ~**±0.002 Å** of pure arithmetic noise. The charter's targets stay
checkable (< 3.00 is 0.21 Å away, < 2.50 is 0.71 Å). **Unpaired cross-job chain claims below ~0.03 Å
are not resolvable**, and "the same cloud value" does not license a comparison — it must be the same
float64 bits.

---

## 2. Where the RMSD is actually lost

[PENDING — the ladder with all order-statistic pricing and stratification.]

---

## 3. The readout is a hull projection, and that closes a family

[PENDING — S32-L5, L(Q1)–L(Q3): gain, the sufficient statistic, the hull floor, monotonicity.]

---

## 4. In-band skill is not zero — the per-target SIGN is missing

[PENDING — S32-L7, L(D2), and whether anything native-free reads the bit.]

---

## 5. The projection taxes direction, not distance

[PENDING — S32-L9, R-13's off-manifold law, the dilation refutations.]

---

## 6. The quantum question, answered

[PENDING — charter §14, the five conditions, P1 and P2, and why chain length is the binding one.]

---

## 7. Longer proteins

[PENDING — lane L, and why the requirement is derived rather than suggested.]

---

## 8. What was falsified, including by its own author

[PENDING — the registered falsifiers that fired, and the corrections.]

---

## 9. Is it possible to lower the RMSD?

[PENDING — the charter's closing question, answered directly.]

---

## 10. The next bottleneck

[PENDING]

---

## Appendix A — every claim withdrawn this sprint

[PENDING]

## Appendix B — multiplicity and the search that was run

[PENDING]
