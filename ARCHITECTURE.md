# ARCHITECTURE — FROZEN 2026-09-08 (Sprint 25)

Authoritative specification of the production pipeline. Frozen at commit `a15406c`.
Anything not described here is research code and lives outside the production path.

**Primary endpoint:** mean full-chain Cα-RMSD over the 126 cluster-disjoint development targets.
**The 60-target benchmark is sealed** and was not used to select any component, parameter or threshold.

---

## 0. THE TWO NUMBERS, AND WHY THERE ARE TWO

| | value | what it is |
|---|---|---|
| **PRODUCTION RESULT** | **3.2148 Å** | `rmsd_arm` — the ideal-geometry chain the system emits. A real, renderable, chemically valid protein structure. |
| intermediate | 3.0483 Å | `rmsd_avg` — the raw coordinate average. **Not a protein structure.** |

The coordinate average is **22.3% contracted**: mean virtual Cα–Cα bond **2.9614 Å** against a native
**3.8122 Å**, worst single bond **0.649 Å** — shorter than a covalent C–C bond — with 80 of 126 targets
under 3.4 Å. `core/bench.py:682` labels this arm `"raw average (illegal)"`. It cannot be written as a
valid PDB, so it cannot round-trip through the structure export, and it is **never** the system's
result.

Projecting it onto the ideal-geometry manifold costs **+0.1664 Å** (median +0.0977, IQR
[+0.005, +0.338]; the chain is *better* on 16 of 126). That cost correlates with each target's own
contraction at **ρ = +0.5179** — on the least-contracted quartile the projection is essentially free.
**It is the price of undoing the averaging's contraction, not a toll for becoming a structure.**

Both numbers were independently re-derived at n=126 through a different RMSD implementation than the
one that wrote them, with **zero disagreement**, and the production arm round-trips from its exported
files at worst 2.641e-04 Å.

> **Point-cloud and built-chain RMSD are different endpoints and are never compared.** Every figure in
> this project states its basis.

---

## 1. DATA FLOW

```
sequence (9–16 aa)
   │
   ├─► ESM-2 embedding ──► distogram model (leave-fold-out, 5 folds)
   │                          │
   │                          ▼
   │                    per-pair distance posterior p(d) over 17 bins
   │                          │
   ├─► BLOSUM window retrieval over the fragment/peptide corpus
   │        │
   │        ▼
   │   K = 500 candidate windows (real protein geometry)
   │        │
   │        ▼
   │   score:  s(w) = Σ_pairs w_p · risk_p(d_p(w))          ← Bayes risk against the posterior
   │        │
   │        ▼
   │   top-128 by score  ──► E = zrank(scores)  ──► H = diag(E)
   │                                                   │
   │                                                   ▼
   │                                    CVaR-VQE  (StatevectorCircuit, 7 qubits,
   │                                     3 layers, 21 params, exact statevector,
   │                                     Adam on exact parameter-shift gradient)
   │                                                   │
   │                                                   ▼
   │                                        p_θ  ──►  CVaR_α tail  ──► selected set
   │        ┌──────────────────────────────────────────┘
   │        ▼
   │   uniform coordinate average of the retained set  ──►  point cloud  (3.0483 Å, NOT a structure)
   │        │
   │        ▼
   │   ideal-geometry projection (Ramachandran-penalised, multi-start)
   │        │
   │        ▼
   └─►  FINAL STRUCTURE — Cα trace of an ideal-geometry backbone   **3.2148 Å**
```

Legacy and AMBER are genuine, independently evaluable Hamiltonians available at the scoring stage.
**Neither is in the production score** — §5 explains why, and it is a measured result, not an omission.

---

## 2. COMPONENTS

### 2.1 Retrieval — `core/peptide_db.py`, `core/data.py`
BLOSUM-scored window search over two disjoint banks: 787 peptide chains (8–26 aa) and 6,003 fragments
(9–20 aa, 1,001 deposits). K = 500 candidates per target. Candidates are **real protein geometry**, so
every candidate is already a valid chain.

**Known defect, declared:** `core/data.py:188` normalises Needleman–Wunsch identity by the *longer*
sequence — its own docstring calls this "leaky at the member level". Consequence: an 11-mer sitting
verbatim inside a 21-mer scores 0.52, falls below the 0.60 clustering threshold, and lands in a
different fold. **4 of 126 dev targets and 2 of 60 benchmark targets carry a verbatim self-copy in
their own fold model's training set and retrieval library.** On the four dev cases the self-window is
BLOSUM rank #1 but is *worse* than the pool's own best on 3 of 4. The caveat attaches to every figure.

### 2.2 Structural prior — `core/predict.py`
ESM-2 sequence embedding → per-pair distance posterior over **17 irregular bins**, centres 4.0 … 25.0 Å
with 0.5 Å spacing at the short end and 3.5 Å at the long end. Leave-fold-out: five models, each
trained without its own fold.

**Measured properties (Sprint 25):** the posterior is **over-confident by roughly 2×** — z_sd 1.66–2.09
across robust estimators against a nominal 1.0, 90% coverage 0.62. Nearly centred. 24.1% of pairs carry
a genuinely multimodal posterior. The defect lives in the **tails**, not the core.

**Note for anyone reading the fold models' independence:** each of the five saw ~100 of the other 125
dev natives. **Per-target dev results are therefore not independent and iid confidence intervals are
anticonservative. Fold-clustered intervals are the only honest ones on this instrument.**

### 2.3 Score — `core/predict.py:420`
```
risk_p(t) = w_p · Σ_c prob[p,c] · |t − CENTRES[c]|          w_p = shell / (sd_p + 0.5)^g
s(w)      = Σ_p risk_p( d_p(w) )
```
`_score_weights()` returns `ones, 1.0` (the weights file is intentionally absent), so **shell ≡ 1 and
g = 1**: the only per-pair weight is a pure function of the posterior's width.

This is an **L1 Bayes risk**, so its per-pair minimiser is the posterior **median** — which, for a
discrete distribution, always sits on an atom. **The score's effective per-pair target therefore takes
only 17 distinct values**, largest adjacent gap 4.0 Å, admitting up to 2.0 Å of pure location error
before any estimation error. Removing that exactly is worth 0.003 Å (§5).

### 2.4 Selector — `core/quantum.py`, `core/pipeline.py:853`
`StatevectorCircuit`, **7 qubits, 3 layers, 21 parameters**, real-amplitude RY/CNOT, exact dense
statevector (verified to 5.6e-17 against an independently written simulator). Adam, 50 iterations, on
the **exact parameter-shift gradient** (cos 1.000000000 against finite differences, rel err 4.6e-10).

Objective: `F(p) = E_p[E] − T·H(p)` with `E = zrank(scores of the top-128)`.
CVaR at level α via `cvar_from_probs`. Deployed table, `core/pipeline.py:118`:
```
VQE_LFO = {0:(1.0,0.3), 1:(0.25,0.3), 2:(1.0,0.3)... }   → T = 0.3 everywhere, α = 1.0 on 3 of 5 folds
```
It is a **leave-fold-out** table: fold *f*'s cell was chosen on the other four.

**`MPSAnsatz` is NOT instantiated by the production pipeline.** It is the generation lane's ansatz.
The MPS is exact (χ = 2^layers by construction, no SVD, no truncation) but it is not the deployed
selector.

### 2.5 Readout — `core/pipeline.py`
Superpose the retained set on its medoid, take the **uniform arithmetic mean**. Uniform weights are
optimal, not merely conventional: a leakage-oracle grid search over rank-power and distance-to-medoid
weightings selects the uniform point, and removing the four most geometrically deviant members costs
+0.142 Å while removing four at random costs nothing.

### 2.6 Projection — `core/project.py`
Nearest ideal-geometry chain with a Ramachandran penalty, multi-start. Emits φ/ψ and the Cα trace of a
real backbone. Bond length is a **constant 3.8040 Å** — idealised geometry, so it cannot represent a
cis-peptide (cis-proline is ~2.9 Å). Stated as a modelling limitation, not a defect.

---

## 3. WHAT THE QUANTUM COMPONENT CONTRIBUTES

**It is genuine.** Exact statevector, exact parameter-shift gradient, a real Hamiltonian-derived
objective, a real CVaR. All re-verified independently in Sprint 25.

**Its selection is provably classical.** The realised CVaR tail's support is always a **subset of an
initial prefix of the energy order**, and equals that prefix exactly when every prefix state carries
positive probability. `p_θ` can *delete* a member of the classical top-*m*; it can never *add* one from
outside it. Verified on 2,592 adversarial cells with exact zeros: 0 subset violations, 972/972 equality
on full support.

**The optimiser genuinely trains** — it beats best-of-200 from the untrained circuit at every
temperature, closing 78–89% of the available free-energy gap.

**And the endpoint cannot see any of it.** RMSD tracks the *entropy of the readout weights*
(ρ = −0.7423), not α (+0.2700) or T (−0.0234). Fit that curve on the **nine arms the circuit is not
in** and the nine circuit arms land on it at +0.0090 Å against a residual sd of 0.0268. The trained
state sits **0.902 nats and 45% of its mass** from its own analytic Gibbs optimum `exp(−E/T)/Z`, and
the endpoint difference is 0.24× MDE — **NULL**.

> **The readout is insensitive to a distributional difference of nearly half the mass.**
> The binding constraint is readout slack. It would bind identically on a landscape where the Gibbs
> state were expensive, so "the target is classical and cheap" is *not* the mechanism.

**And the Hamiltonian barely changes between targets.** Because `top = argsort(sc)`, `E` is the
standardised rank ladder up to tie-averaging — worst deviation 1.18% of range across targets. There are
effectively **two trained states in the whole deployment, not 126**, and all per-target information
enters through which candidate occupies which rank. **This is why deeper ansätze and larger χ order
nothing: extra expressivity has nothing to be expressive about.**

---

## 4. NORMALISATION BETWEEN HAMILTONIANS

`zrank(x) = (rankdata(x) − mean)/sd`; with no ties `mean = (K+1)/2` and `sd = √((K²−1)/12)` exactly, so
every channel lands on the identical marginal. Combined channels use
`E_S = zrank(Σ_{c∈S} zrank(x_c))` — **the outer re-standardisation is load-bearing**, because CVaR
trades energy against `T·H` and without it the effective temperature would differ silently between
configurations (worth up to 0.14 Å).

**Raw moment standardisation is not usable here.** AMBER on unrelaxed retrieval windows is
clash-dominated — 57.5% of candidates above 1e4 kcal/mol, worst 7.1e18 — which puts 99.4% of a pool
inside |z| < 0.1. Worse, `zmoment` is **not monotone in float64** when the sd is set by a 1e28 outlier:
it breaks AMBER's own `argsort` on 40 of 126 targets, and `argsort` then silently reads the BLOSUM
retrieval order, which is not neutral.

---

## 5. KNOWN LIMITATIONS

1. **<3.0 Å was not reached.** The production result is 3.2148 Å built chain / 3.0483 Å point cloud.
2. **The binding constraint is the distance prior's accuracy, and it is not reachable by re-reading.**
   An oracle ladder that interpolates the posterior toward the truth moves the endpoint at −2.15 Å per
   unit at the origin, and γ = 0.0225 would reach 3.0 Å — **but that ladder's γ is redeemable only
   along the native's own direction (cos = 1 by construction).** A real operator travelling 25% of the
   way at cos = 0.50 is worth +0.024 Å. Never quote −2.1496 × γ for an achievable operator without the
   cosine beside it.
3. **Every re-reading of the existing posterior is measured out.** Width calibration, location shifts
   (global and separation-graded), mode-seeking, alternative risk functionals, metric projection and
   de-quantisation: for three of four families the best parameter chosen with *complete leakage* is the
   shipped identity. Across 51 arms, `corr(progress-toward-truth, endpoint) = +0.054`.
4. **Both physics energies are worse than noise as rankers on this pool** (§L16). Legacy +0.330 and
   AMBER +0.455 against a random 75-subset, 5/5 folds. Permuting a physics channel while preserving its
   marginal *improves* the endpoint.
5. **Candidate generation is closed** on five instruments (Sprint 24), including a learned residual
   generator closed by an oracle upper bound and a from-scratch generator closed by a corpus census.
6. **iid CIs across this project are anticonservative** — see §2.2.
7. **The 4/126 dev and 2/60 benchmark self-copy leak** — see §2.1. The benchmark remains sealed; the
   2/60 is declared and was never quantified, because quantifying it requires benchmark RMSD.

---

## 6. REPRODUCING

```
commit           a15406c
instrument       s12/instrument.py — K=500, M=75, 126 targets in pinned pdb-sorted order, 5 pinned folds
corpus hash      29e3b67e8ca0c03d
seeds            s15/seed.py stable_rng(experiment, key)
statistics       s24/stats_lib.py — MDE = 2.8016 × SE per comparison, fold-clustered CI beside iid,
                 best_of_k_within for K-settings grids, split_half_transfer, achieved()
artefacts        ST.save_atomic(..., module_file=__file__) — module sha256, git commit, dirty flag
```

**Release gates for any generated result set:** mean RMSD ≥ pool_best 1.7108 Å; `corr(per-target RMSD,
target difficulty)` and per-target sd above floors; a **mandatory provenance REMARK in every PDB
header** with the build refusing any file lacking one; and test fixtures writing to a temp directory,
never the real output tree.
