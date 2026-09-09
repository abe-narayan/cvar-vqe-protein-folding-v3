"""The redesigned folding pipeline: prior -> warm-started CVaR-VQE -> refine -> select.

Pipeline
::

    sequence
       |
       +-- peptide_db.holdout ------------------------------------------+
       |      (target and everything >0.6 identity to it removed)       |
       v                                                                |
    per-residue torsion library (torsion_lib2)  <---------- context stats
       |                                                                |
       +--> torsion MRF prior  --------+                                |
       |        (h_i, J_i,i+1)         |                                |
       +--> ESM-2 + sequence distogram-+                                |
                |                      |                                |
                v                      v                                |
          FoldObjective         warm-start angles ----+                 |
                |                                     |                 |
                v                                     v                 |
        CVaR-VQE  (exact-sampled RY/CNOT ansatz, analytic CVaR gradient) |
                |                                                       |
                +--> diverse reservoir of candidates                    |
                             |                                          |
                             v                                          |
                  continuous torsion refinement                         |
                             |                                          |
                             v                                          |
              basin selection (population + score)                      |
                             |                                          |
                             v                                          |
              optional physics rescoring: Model A (energy_terms)        |
                                          Model B (Amber ff14SB + GB) --+

What changed and why
*Warm start.* The RY angles of the one-layer ansatz set the per-qubit marginals exactly,
and the marginals of a torsion-state prior are exactly what we know before searching. The
inversion is closed form: for the CNOT chain, ``P(b_q = 1) = (1 - prod_{t<=q}(1 - 2 s_t))
/ 2``, so the required Bernoulli parameters follow from a running product. The previous
driver initialised at pi/2 plus noise -- the uniform distribution -- and spent its budget
rediscovering the Ramachandran statistics of a 780-peptide database.

*Real CVaR gradient.* `qansatz.cvar_gradient` differentiates the CVaR objective through
the score function instead of pulling the angles toward the elite mean bit pattern.

*Adaptive alpha.* Annealed from 0.5 to 0.05 across a restart: broad while the landscape is
still being mapped, narrow once it is.

*Selection by basin, not by score.* The scorer's own ranking of its top candidates is
weak; the population of a basin under a decent scorer is a much better signal, and it is
what guards the worst case.
"""
import os
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import distogram as dgm
import objective as obj
import peptide_db as pdb
import priors
import protein_geometry as geo
import qansatz
import torsion_lib2 as tl2


# --------------------------------------------------------------------- warm start
def marginals_to_angles(bit_marginals: np.ndarray, ring: bool = False) -> np.ndarray:
    """RY angles whose one-layer ansatz reproduces the given per-qubit marginals.

    ``b_q`` is the prefix XOR of independent Bernoulli(s_t), so
    ``1 - 2 P(b_q = 1) = prod_{t <= q} (1 - 2 s_t)``. Dividing consecutive running
    products recovers each ``s_q`` exactly. The ring closure is ignored here: it couples
    the first and last qubits, and matching marginals up to that one coupling is enough
    for an initialisation.
    """
    m = np.clip(np.asarray(bit_marginals, float), 0.02, 0.98)
    A = 1.0 - 2.0 * m
    prev = np.concatenate([[1.0], A[:-1]])
    ratio = np.clip(A / prev, -0.999, 0.999)
    s = np.clip((1.0 - ratio) / 2.0, 1e-3, 1 - 1e-3)
    return 2.0 * np.arcsin(np.sqrt(s))


def state_marginals_to_bits(pi: np.ndarray, bits_per_residue: int) -> np.ndarray:
    """``pi`` (n, k) state probabilities -> (n * bits_per_residue,) qubit marginals."""
    n, k = pi.shape
    codes = ((np.arange(k)[:, None] >> np.arange(bits_per_residue - 1, -1, -1)) & 1)
    return (pi @ codes).reshape(-1)


def make_ansatz(kind: str, n_qubits: int, layers: int, rep, mrf,
                warm_start: bool):
    """``(ansatz, initial angles)``.

    ``chain``     the shipped one-layer RY + CNOT chain + ring, sampled in closed form.
    ``chain_ry``  the same chain with one trailing RY layer, which is what makes the
                  prior's marginals reachable at all (see `qansatz.MPSAnsatz`).
    ``deep``      ``layers`` chain blocks plus a trailing RY layer, exact MPS.
    ``product``   RY only, no entangler -- the control for whether entanglement helps.
    """
    m = (state_marginals_to_bits(mrf.marginals(), rep.bits_per_residue)
         if warm_start else None)
    warm = marginals_to_angles(m) if warm_start else np.full(n_qubits, np.pi / 2.0)
    if kind == "chain":
        return qansatz.OneLayerAnsatz(n_qubits, ring=True), warm
    if kind == "product":
        return (qansatz.MPSAnsatz(n_qubits, 1, final_ry=False, entangler="none"),
                warm)
    nblocks = 1 if kind == "chain_ry" else max(1, layers)
    a = qansatz.MPSAnsatz(n_qubits, nblocks, final_ry=True)
    # Warm start goes in the LAST rotation layer with the entangling blocks at identity:
    # RY(0) is the identity and a CNOT on |0...0> does nothing, so the prepared state is
    # exactly the product distribution with the prior's marginals. Putting it in the first
    # layer instead would push it through the CNOT chain and lose it.
    base = np.concatenate([np.zeros(n_qubits * nblocks), warm])
    return a, base


# --------------------------------------------------------------------- reservoir
class Reservoir:
    """Score-ordered, diversity-preserving candidate store."""

    def __init__(self, radius: float = 1.2, capacity: int = 64):
        self.radius = radius
        self.capacity = capacity
        self.scores: List[float] = []
        self.cas: List[np.ndarray] = []
        self.states: List[np.ndarray] = []

    def offer_batch(self, scores: np.ndarray, cas: np.ndarray,
                    states: np.ndarray) -> None:
        for s, c, st in zip(scores, cas, states):
            self._offer(float(s), c, st)

    def _offer(self, score: float, ca: np.ndarray, st: np.ndarray) -> None:
        if not self.cas:
            self.scores.append(score)
            self.cas.append(ca)
            self.states.append(st)
            return
        d = geo.ca_rmsd_batch(np.stack(self.cas), ca)
        j = int(np.argmin(d))
        if d[j] > self.radius:
            self.scores.append(score)
            self.cas.append(ca)
            self.states.append(st)
            if len(self.scores) > self.capacity:
                worst = int(np.argmax(self.scores))
                for lst in (self.scores, self.cas, self.states):
                    lst.pop(worst)
        elif score < self.scores[j]:
            self.scores[j], self.cas[j], self.states[j] = score, ca, st

    def items(self):
        return list(zip(self.scores, self.cas, self.states))


# --------------------------------------------------------------------- refinement
def refine(objective, states: np.ndarray, sigma0: float = 12.0, steps: int = 40,
           pop: int = 24, seed: int = 0, restarts: int = 2
           ) -> Tuple[np.ndarray, np.ndarray, float]:
    """Continuous local relaxation of (phi, psi) away from the discrete grid.

    A (1 + pop) evolution strategy in torsion space with a geometrically shrinking step,
    which reuses the batched builder that the search already runs at ~4,000 structures per
    second. The discrete grid costs ~0.5-1.5 A of representation floor; this is what
    recovers the part of that gap the objective can actually see.
    """
    rows = np.arange(objective.n)
    phi0 = objective.rep._phi[rows, states].copy()
    psi0 = objective.rep._psi[rows, states].copy()
    e0, ca0 = objective.score_angles(phi0[None], psi0[None])
    gbest, gphi, gpsi, gca = float(e0[0]), phi0, psi0, ca0[0]
    # Several short anneals from the same discrete start beat one long one: the surface off
    # the grid is rough enough that a single (1+lambda) run stalls at the first basin it
    # shrinks into.
    for r in range(max(1, restarts)):
        rng = np.random.default_rng(seed * 977 + r)
        phi, psi = phi0.copy(), psi0.copy()
        best, best_ca = float(e0[0]), ca0[0]
        sigma = np.radians(sigma0)
        for t in range(steps):
            cand_phi = phi[None] + rng.normal(0, sigma, (pop, objective.n))
            cand_psi = psi[None] + rng.normal(0, sigma, (pop, objective.n))
            e, cas = objective.score_angles(cand_phi, cand_psi)
            k = int(np.argmin(e))
            if e[k] < best:
                best, phi, psi, best_ca = float(e[k]), cand_phi[k], cand_psi[k], cas[k]
            else:
                sigma *= 0.88
            if sigma < np.radians(0.4):
                break
        if best < gbest:
            gbest, gphi, gpsi, gca = best, phi, psi, best_ca
    return gphi, gpsi, gbest, gca


# --------------------------------------------------------------------- selection
def snap_to_states(rep, ca: np.ndarray, seed: int = 0, restarts: int = 6):
    """Torsion states whose CA trace best matches `ca`, and that trace."""
    import floor
    rng = np.random.default_rng(seed)
    best, best_v = None, float("inf")
    for r in range(restarts):
        start = (np.zeros(rep.n_residues, int) if r == 0
                 else rng.integers(0, rep.n_states, rep.n_residues))
        st, v = floor.descend(rep, np.asarray(ca, float), start)
        if v < best_v:
            best, best_v = st, v
    rows = np.arange(rep.n_residues)
    out = geo.build_backbone_batch(rep._phi[rows[None], best[None]],
                                   rep._psi[rows[None], best[None]])["CA"][0]
    return best, out


def basin_average(cas: Sequence[np.ndarray], weights: Sequence[float],
                  members: Sequence[int], centre: int) -> np.ndarray:
    """Boltzmann-weighted mean of a basin's members, superposed on its medoid.

    Each member carries independent error; averaging them cancels the part of it that is
    random rather than systematic. The result is not a lattice structure -- it does not
    correspond to any torsion assignment -- so it is only ever the *reported* answer, never
    a search candidate.
    """
    ref = np.asarray(cas[centre], float)
    w = np.asarray([weights[i] for i in members], float)
    w = w / max(w.sum(), 1e-12)
    stack = np.stack([geo.kabsch_superpose(np.asarray(cas[i], float), ref)
                      for i in members])
    return (stack * w[:, None, None]).sum(0)


def basin_select(scores: Sequence[float], cas: Sequence[np.ndarray],
                 radius: float = 2.0, kT: float = 1.0) -> Tuple[int, Dict]:
    """Medoid of the basin with the most Boltzmann-weighted population.

    Population, not best score: on this objective the top candidates cannot be ranked
    against each other reliably, but a basin that many independent restarts fall into is
    much more often the right one. `kT` is in units of the objective, whose scale is set
    by the distogram term (mean per-pair distance error in Angstroms).
    """
    C = np.stack(cas)
    s = np.asarray(scores, float)
    D = np.stack([geo.ca_rmsd_batch(C, c) for c in C])
    w = np.exp(-(s - s.min()) / kT)
    mass = (D < radius) @ w
    members = np.where(D[int(np.argmax(mass))] < radius)[0]
    sub = D[np.ix_(members, members)]
    pick = int(members[np.argmin((sub * w[members][None, :]).sum(1))])
    return pick, {"n_basins": int(len(set(map(tuple, (D < radius).astype(int))))),
                  "basin_size": int(len(members)),
                  "basin_mass": float(mass.max() / w.sum()),
                  "members": members.tolist(), "weights": w.tolist()}


# --------------------------------------------------------------------- driver
def build_distogram(sequence: str, rep, use_esm: bool = True, fragments: bool = True,
                    models: str = "auto", n_ref: int = 512, seed: int = 0):
    """The distance prior for a target: one model, or both combined and calibrated.

    ``models="both"`` averages the per-pair MLP and the pair-tensor network after
    standardising each against a fixed sample of random structures from this target's own
    state library. The sample is drawn from the representation, never from the native.
    """
    if models == "auto":
        # "both" is better (in-band +0.454 against +0.451 and +0.397 alone), but only when
        # the pair-network fold models actually exist. Training one on demand inside a
        # target build is a 25-minute stall, so availability is checked rather than assumed.
        import glob
        import pairnet as _pn
        models = "both" if glob.glob(os.path.join(_pn.MODEL_DIR, "fold*.pt")) else "mlp"
    parts = []
    if models in ("mlp", "both"):
        parts.append(dgm.Distogram.for_target(sequence, use_esm=use_esm,
                                              fragments=fragments))
    if models in ("pairnet", "both"):
        try:
            import pairnet
            parts.append(pairnet.distogram_for(sequence))
        except Exception:
            pass
    if len(parts) == 1:
        return parts[0]
    rng = np.random.default_rng(seed)
    S = rng.integers(0, rep.n_states, size=(n_ref, rep.n_residues))
    rows = np.arange(rep.n_residues)
    ref = geo.build_backbone_batch(rep._phi[rows[None], S], rep._psi[rows[None], S])["CA"]
    return dgm.CombinedDistogram(parts, ref_ca=ref)


def build_target(sequence: str, k: int = 8, use_esm: bool = True,
                 exclude_self: bool = True, w_mrf: float = 0.35,
                 w_clash: float = 1.0, fragments: bool = True,
                 models: str = "auto", w_catrace: float = 0.0,
                 w_legacy: float = 0.0,
                 legacy_weights: Optional[Dict[str, float]] = None,
                 legacy_use: Optional[Sequence[str]] = None) -> Dict:
    """Everything that depends only on the sequence: library, priors, objective."""
    ex = sequence if exclude_self else ""
    entries = pdb.holdout(sequence) if exclude_self else list(pdb.load())
    tab = tl2.library_for(sequence, k, ex)
    rep = tl2.PerResidueTorsion(sequence, tab, chi_bits=False)
    dist = build_distogram(sequence, rep, use_esm=use_esm, fragments=fragments,
                           models=models)
    mrf = priors.TorsionMRF.fit(sequence, tab, entries)
    import catrace
    catr = catrace.CATracePrior(exclude_seq=ex)
    leg = None
    if w_legacy:
        import legacy_field
        leg = legacy_field.LegacyField(sequence, rep, weights=legacy_weights,
                                       use=legacy_use, dist=dist)
    o = obj.FoldObjective(sequence, rep, dist=dist, mrf=mrf, catr=catr,
                          w_mrf=w_mrf, w_clash=w_clash, w_catrace=w_catrace,
                          legacy=leg, w_legacy=w_legacy)
    return {"rep": rep, "dist": dist, "mrf": mrf, "catrace": catr, "legacy": leg,
            "objective": o, "table": tab}


def seed_states(target: Dict, restarts: int = 4, seed: int = 0) -> np.ndarray:
    """Torsion states whose CA trace best matches the distogram's own embedded structure.

    Distance geometry on the *predicted* matrix gives a coordinate trace; snapping that
    trace onto the state library is the discrete assignment the search should start from.
    Nothing here touches the native -- the target's structure is not read at any point --
    so this is a prior-driven initialisation, not leakage. `floor.descend` is reused
    verbatim: it minimises CA-RMSD of a state assignment against a given trace, and here
    the given trace is the prediction rather than the answer.
    """
    import floor
    rep = target["rep"]
    X = target["dist"].realize(seed=seed, restarts=restarts)
    rng = np.random.default_rng(seed)
    best, best_v = None, float("inf")
    for r in range(6):
        start = (np.zeros(rep.n_residues, int) if r == 0
                 else rng.integers(0, rep.n_states, rep.n_residues))
        st, v = floor.descend(rep, X, start)
        if v < best_v:
            best, best_v = st, v
    return best


def run_vqe(objective, ansatz, theta0: np.ndarray, iters: int, shots: int,
            rng, lr: float = 0.12, a0: float = 0.5, a1: float = 0.05,
            reservoir: Optional[Reservoir] = None, elite_frac: float = 0.1,
            random_frac: float = 0.0,
            update: str = "cvar_grad", trace: Optional[List] = None) -> Dict:
    """One CVaR-VQE restart. Returns the final angles and the best structure seen."""
    theta = np.array(theta0, float)
    adam = qansatz.Adam(len(theta), lr=lr)
    rep = objective.rep
    w, n_res, k = rep.bits_per_residue, rep.n_residues, rep.n_states
    powers = (1 << np.arange(w - 1, -1, -1)).astype(np.int64)
    best_e, best_ca, best_st = float("inf"), None, None
    for it in range(iters):
        alpha = qansatz.alpha_schedule(it / max(1, iters - 1), a0, a1)
        bits = ansatz.sample(theta, shots, rng)
        states = (bits[:, :w * n_res].reshape(shots, n_res, w) @ powers) % k
        e, ca = objective(states)
        m = int(np.argmin(e))
        if e[m] < best_e:
            best_e, best_ca, best_st = float(e[m]), ca[m], states[m]
        if reservoir is not None:
            sel = np.argsort(e)[:max(1, int(elite_frac * shots))]
            if random_frac > 0:
                # Admission by SCORE is admission by a weak ranker. Instrumentation shows
                # the best structure the objective ever scores averages 1.06 A while the
                # best one that survives the elite filter averages 1.23 A -- the filter
                # itself throws away 0.17 A using a score whose in-band Spearman is 0.45.
                # Offering a score-blind random sample as well lets the diversity radius,
                # rather than the score, decide what the pool covers.
                extra = rng.choice(shots, size=max(1, int(random_frac * shots)),
                                   replace=False)
                sel = np.unique(np.concatenate([sel, extra]))
            reservoir.offer_batch(e[sel], ca[sel], states[sel])
        if update == "cvar_grad":
            grad, val = qansatz.cvar_gradient(ansatz, theta, bits, e, alpha)
            theta = adam.step(theta, grad)
        else:
            # The shipped driver's update, kept so the comparison is measurable: pull the
            # angles toward the mean bit pattern of the CVaR tail. Not a gradient of any
            # objective -- it is the cross-entropy method.
            val = qansatz.cvar(e, alpha)[0]
            tail = np.argsort(e)[:max(1, int(np.ceil(alpha * shots)))]
            frac = np.clip(bits[tail].mean(0), 1e-6, 1 - 1e-6)
            tgt = 2.0 * np.arcsin(np.sqrt(frac))
            step_lr = lr * 2.9 * (1.0 - it / max(1, iters - 1))
            theta = theta + step_lr * (tgt - theta) + rng.normal(0, 0.05, theta.shape)
        if trace is not None:
            trace.append({"iter": it, "alpha": alpha, "cvar": val,
                          "min": float(e.min()), "mean": float(e.mean()),
                          "best": best_e})
    return {"theta": theta, "best_energy": best_e, "best_ca": best_ca,
            "best_states": best_st}


def fold(sequence: str, k: int = 8, layers: int = 1, restarts: int = 4,
         iters: int = 60, shots: int = 384, seed: int = 0, use_esm: bool = True,
         warm_start: bool = True, refine_top: int = 8, radius: float = 1.2,
         fragments: bool = True, select_radius: float = 2.0, target: Optional[Dict] = None,
         exclude_self: bool = True, lr: float = 0.12,
         alpha0: float = 0.5, alpha1: float = 0.05, update: str = "cvar_grad",
         refine_kw: Optional[Dict] = None, random_frac: float = 0.0,
         reservoir_capacity: int = 64,
         select: str = "basin", ansatz_kind: str = "chain",
         seed_frac: float = 0.6, seed_restarts: int = 1, use_seed: bool = True,
         verbose: bool = False) -> Dict:
    """Fold `sequence`. `exclude_self=True` is the benchmark setting: the target's own
    structure and every homolog above 0.6 identity are removed from every fitted
    component."""
    t_start = time.perf_counter()
    T = target if target is not None else build_target(
        sequence, k=k, use_esm=use_esm, exclude_self=exclude_self,
        fragments=fragments)
    rep, o, mrf = T["rep"], T["objective"], T["mrf"]
    o.n_calls = o.n_structures = 0      # per-fold accounting, not cumulative
    n_qubits = rep.bits_per_residue * rep.n_residues

    # Distance-geometry seed: embed the predicted distogram, snap the trace onto the state
    # library, and mix that assignment into the warm-start marginals. The VQE still does
    # the search -- this only says where to start looking, and it is derived entirely from
    # the prediction, never from the target's structure.
    seed_st = None
    if use_seed:
        try:
            seed_st = seed_states(T, seed=seed)
        except Exception:
            seed_st = None

    ansatz, base = make_ansatz(ansatz_kind, n_qubits, layers, rep, mrf, warm_start)
    seeded_base = base
    if seed_st is not None and seed_frac > 0 and warm_start:
        pi = mrf.marginals()
        onehot = np.zeros_like(pi)
        onehot[np.arange(len(seed_st)), seed_st] = 1.0
        blended = (1.0 - seed_frac) * pi + seed_frac * onehot
        seeded_base = make_ansatz(
            ansatz_kind, n_qubits, layers, rep,
            priors.TorsionMRF(-np.log(np.clip(blended, 1e-6, None)), mrf.J),
            True)[1]

    res = Reservoir(radius=radius, capacity=reservoir_capacity)
    if seed_st is not None:
        e_seed, ca_seed = o(seed_st[None])
        res.offer_batch(e_seed, ca_seed, seed_st[None])
    traces, runs = [], []
    for r in range(restarts):
        rng = np.random.default_rng([seed, r])
        # The first `seed_restarts` restarts start from the distance-geometry seed's
        # basin; the rest start from the sequence prior. Seeding every restart measurably
        # narrows the search (it wins on some targets and loses on more), so the two
        # initialisations are run side by side inside one fold instead of chosen.
        start = seeded_base if r < seed_restarts else base
        theta0 = start + rng.normal(0, 0.35, size=ansatz.n_params())
        tr: List = []
        run = run_vqe(o, ansatz, theta0, iters, shots, rng, lr=lr,
                      a0=alpha0, a1=alpha1, reservoir=res, update=update,
                      random_frac=random_frac, trace=tr)
        runs.append(run)
        traces.append(tr)
        if verbose:
            print(f"    restart {r}: best {run['best_energy']:.4f}, "
                  f"reservoir {len(res.scores)}", flush=True)

    items = res.items()
    scores = [s for s, _, _ in items]
    cas = [c for _, c, _ in items]
    states = [st for _, _, st in items]
    # Torsion angles of every pool member, kept so a *resolution-independent* rescore can
    # use terms that need more than the CA trace. A refined candidate no longer sits on
    # its state library, so its states are not a description of it; its angles are.
    _rows = np.arange(rep.n_residues)
    _S = np.asarray(states, int)
    phis = list(rep._phi[_rows[None, :], _S]) if len(states) else []
    psis = list(rep._psi[_rows[None, :], _S]) if len(states) else []

    # Refinement. Two corrections over the previous version, both from measurement:
    #
    # 1. The old code compared a refined `dist + clash` value against an unrefined
    #    `dist + MRF + clash` one, so a refinement was accepted partly because the MRF term
    #    silently vanished. That mismatch was acting as an accidental conservative accept
    #    filter -- putting both on one basis measures 0.14 A WORSE. Rather than keep an
    #    accident, the refined structure is now ADDED to the pool alongside the unrefined
    #    one and everything is re-scored on the same basis, so selection sees both and no
    #    hidden filter decides for it.
    # 2. `refine2.refine_pool` relaxes the whole reservoir with an exact analytic torsion
    #    gradient at 0.33 s/structure -- 2.7x cheaper than the (1+lambda) search it
    #    replaces, and it reaches a lower objective. Accuracy is the same, because the
    #    objective, not the operator, is the constraint: across nine refinement settings
    #    the within-target Spearman between objective reached and CA-RMSD is +0.014.
    refined = None
    if refine_top:
        order = list(np.argsort(scores)[:refine_top])
        new_cas, new_states, new_scores = [], [], []
        new_phi, new_psi = [], []
        try:
            import refine2
            phi_r, psi_r, e_r, ca_r = refine2.refine_pool(
                o, [states[i] for i in order], seed=seed, **(refine_kw or {}))
            for m, idx in enumerate(order):
                new_cas.append(ca_r[m])
                new_states.append(states[idx])
                new_scores.append(float(e_r[m]))
                new_phi.append(phi_r[m])
                new_psi.append(psi_r[m])
        except Exception:
            for idx in order:
                phi, psi, e, ca_r = refine(o, states[idx], seed=seed + int(idx))
                new_cas.append(ca_r)
                new_states.append(states[idx])
                new_scores.append(float(e))
                new_phi.append(phi)
                new_psi.append(psi)
        # Re-score refined structures on the same basis the unrefined ones carry, so the
        # two populations are comparable.
        if new_cas:
            base = o.w_mrf * o.mrf.score(np.stack(new_states)) if (o.mrf is not None
                                                                   and o.w_mrf) else 0.0
            cas.extend(new_cas)
            states.extend(new_states)
            scores.extend(list(np.asarray(new_scores) + base))
            phis.extend(new_phi)
            psis.extend(new_psi)
        refined = int(len(order))

    if select == "bag":
        # Rank the pool by agreement across random pair subsets, then take the medoid of
        # the best-ranked basin. Combines the two selection signals that were measured to
        # matter: basin population and robustness of the score to which pairs it uses.
        rk = T["dist"].bagged_rank(np.stack(cas), n_bags=16, seed=seed)
        pick, info = basin_select(list(rk), cas, radius=select_radius, kT=6.0)
        info["mode"] = "bag"
    elif select == "best":
        pick, info = int(np.argmin(scores)), {"mode": "best_score"}
    elif select == "consensus":
        C = np.stack(cas)
        D = np.stack([geo.ca_rmsd_batch(C, c) for c in C])
        pick, info = int(np.argmin(D.mean(1))), {"mode": "consensus"}
    else:
        pick, info = basin_select(scores, cas, radius=select_radius)
        info["mode"] = "basin"
    wall = time.perf_counter() - t_start
    rows = np.arange(rep.n_residues)
    _st = np.asarray(states[pick], int)
    sel_coords = geo.build_backbone_batch(rep._phi[rows[None], _st[None]],
                                          rep._psi[rows[None], _st[None]])
    sel_coords = {k: v[0] for k, v in sel_coords.items()}
    return {"ca": cas[pick], "score": scores[pick], "states": states[pick],
            "coords": sel_coords,
            "dist": T["dist"], "objective": o, "rep": rep,
            "pool_k": [k] * len(cas), "reps": {k: rep},
            "clash_of": (lambda C: o.w_clash * o.clash({"CA": np.asarray(C, float)})),
            "pool_scores": scores, "pool_cas": cas, "pool_states": states,
            "pool_phi": phis, "pool_psi": psis,
            "best_score_index": int(np.argmin(scores)),
            "selection": info, "n_candidates": len(items),
            "n_structures": o.n_structures, "wall": wall,
            "structures_per_s": o.n_structures / max(wall, 1e-9),
            "qubits": n_qubits, "layers": layers, "k": k,
            "traces": traces, "refined": refined,
            "vqe_best_energy": float(min(r["best_energy"] for r in runs))}


# --------------------------------------------------------------------- multi-resolution
def fold_multi(sequence: str, ks: Sequence[int] = (8, 16), seed: int = 0,
               targets: Optional[Dict[int, Dict]] = None,
               select_radius: float = 2.0,
               w_legacy_select: Optional[float] = None, **kw) -> Dict:
    """Run the pipeline at several state resolutions and select across the union.

    The ablation says k=8 and k=16 are not ordered: k=16 gives a better median and more
    sub-2 A results (its representation floor is 0.69 A against 0.93 A) and a slightly
    worse mean, because the larger register is harder for the same shot budget. Running
    both and selecting over the merged reservoir takes the better of the two per target
    without having to pick in advance.

    Candidates are re-scored on a resolution-independent basis before selection -- the
    distance prior plus the clash term, both functions of the CA trace alone -- because
    the MRF term's scale depends on the library and would otherwise bias the merge toward
    whichever k produced it.
    """
    pools_ca: List[np.ndarray] = []
    pools_st: List[np.ndarray] = []
    pools_k: List[int] = []
    pools_phi: List[np.ndarray] = []
    pools_psi: List[np.ndarray] = []
    per_k: Dict[int, Dict] = {}
    ref = None
    for k in ks:
        T = (targets or {}).get(k)
        r = fold(sequence, k=k, seed=seed, target=T, select_radius=select_radius, **kw)
        per_k[k] = r
        pools_ca.extend(r["pool_cas"])
        pools_st.extend(r["pool_states"])
        pools_phi.extend(r.get("pool_phi", []))
        pools_psi.extend(r.get("pool_psi", []))
        pools_k.extend([k] * len(r["pool_cas"]))
        if ref is None:
            ref = r
    o = ref["objective"] if "objective" in ref else None
    dist = per_k[ks[0]]["dist"]
    C = np.stack(pools_ca)
    scores = list(dist.score(C) + ref["clash_of"](C))
    # The Legacy field, if the search used one, is added on the same resolution-independent
    # basis. It is a function of the torsion ANGLES, not of the state library, so it merges
    # across k exactly as the distance term does -- and unlike the distance term it can see
    # backbone hydrogen bonding, which is why the merged pool can contain a structure the
    # distance score alone has no way to prefer.
    #
    # `w_legacy_select` defaults to the weight the SEARCH used. Setting it to 0 separates
    # the two effects cleanly: the pool is whatever generation produced either way, and
    # only the pick changes.
    leg = getattr(o, "legacy", None) if o is not None else None
    wl = (getattr(o, "w_legacy", 0.0) if w_legacy_select is None
          else float(w_legacy_select))
    if leg is not None and wl and len(pools_phi) == len(pools_ca):
        PHI, PSI = np.stack(pools_phi), np.stack(pools_psi)
        coords = geo.build_backbone_batch(PHI, PSI)
        scores = list(np.asarray(scores)
                      + wl * leg.score_from_coords(coords, phi=PHI, psi=PSI))
    pick, info = basin_select(scores, pools_ca, radius=select_radius)
    wall = sum(r["wall"] for r in per_k.values())
    out = dict(per_k[ks[0]])
    out.update({"ca": pools_ca[pick], "score": scores[pick], "states": pools_st[pick],
                "pool_scores": scores, "pool_cas": pools_ca, "pool_states": pools_st,
                "pool_phi": pools_phi, "pool_psi": pools_psi,
                "pool_k": pools_k, "reps": {k: per_k[k]["rep"] for k in ks},
                "n_candidates": len(pools_ca), "selection": info, "wall": wall,
                "n_structures": sum(r["n_structures"] for r in per_k.values()),
                "structures_per_s": sum(r["n_structures"] for r in per_k.values()) / max(wall, 1e-9),
                "ks": list(ks)})
    return out


# --------------------------------------------------------- iterated consensus pass
def consensus_distances(pool_cas: Sequence[np.ndarray], pool_scores: Sequence[float],
                        dist, kT: float = 1.0) -> np.ndarray:
    """Score-weighted mean CA-CA distance over the pool, indexed like ``dist.i/j``.

    Measured over 16 targets, this beats the PRIOR's own predicted distances:
    correlation with the native matrix 0.698 -> 0.731 and mean absolute error
    2.317 -> 2.217 A. Uniform and Boltzmann weighting are within 0.01 of each other and a
    50/50 blend with the prior is slightly worse than the pure ensemble, so the pool is
    not merely echoing the prior -- it holds distance information the prior does not.
    """
    C = np.stack(pool_cas)
    s = np.asarray(pool_scores, float)
    w = np.exp(-(s - s.min()) / max(kT, 1e-9))
    w = w / w.sum()
    d = np.linalg.norm(C[:, dist.i] - C[:, dist.j], axis=-1)
    return (d * w[:, None]).sum(0)


def fold_iterated(sequence: str, ks: Sequence[int] = (8, 16), seed: int = 0,
                  targets: Optional[Dict[int, Dict]] = None, w_cons: float = 0.35,
                  passes: int = 2, kT: float = 1.0, **kw) -> Dict:
    """Search, build a consensus distance matrix from the pool, then search again.

    The second pass adds an L1 restraint toward the first pass's consensus distances. The
    restraint weight is deliberately small: this is also the obvious way to amplify the
    prior's own mistakes, and `w_cons` is fitted on the development peptides, never on a
    benchmark target.
    """
    T = targets or {k: build_target(sequence, k=k) for k in ks}
    dist = T[ks[0]]["dist"]
    out = None
    for it in range(max(1, passes)):
        out = fold_multi(sequence, ks=ks, seed=seed + 977 * it, targets=T, **kw)
        if it == passes - 1:
            break
        cons = consensus_distances(out["pool_cas"], out["pool_scores"], dist, kT=kT)
        for k in ks:
            T[k]["objective"].set_consensus(cons, dist.w, w_cons)
    for k in ks:
        T[k]["objective"].set_consensus(None)
    out["passes"] = passes
    out["w_cons"] = w_cons
    return out
