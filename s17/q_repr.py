"""SPRINT 17 / QUANTUM -- T2, ENSEMBLE REWEIGHTING (sprint section 24) and ALPHA-AS-TEMPERATURE.

PRE-REGISTERED IN `s17/PREREG_quantum.md` (P4).  Fires P4 and supplies the exact half of P6.

THE ONE THING A TRAINED CIRCUIT OWNS THAT NO CLASSICAL SEARCH DOES is a normalised
distribution `q_theta(x)` defined on the WHOLE register -- including configurations it never
drew.  Sprint 16 only ever consumed the VQE's *samples*.  Section 24 asks the different
question: use the VQE to estimate a WEIGHTING over candidates that already exist.

HYPOTHESIS (P4).  `log q_theta` carries ranking information about true RMSD BEYOND the
objective `E` it was trained on -- i.e. the ansatz's bounded correlation structure acts as a
useful regulariser of a noisy objective.

EXPECTED OUTCOME.  No.  A Boltzmann law of `E` has EXACTLY the same ranking as `-E` (a
monotone map), so a "Boltzmann reweight" is not a distinct ranker at all; the only thing that
can differ is the SMOOTHING a restricted variational family imposes.  The classical control
for that is the fully factorised (mean-field) Boltzmann law of `E`, computed exactly here.

CONTROLS, both mandatory.
  * `meanfield` -- the exact zero-correlation limit of the variational family, at matched
    entropy.  This is the correct classical analogue and Sprint 16 never ran it.
  * `-E` itself -- the zero-information-added reference: whatever `q_theta` does, it must
    beat simply ranking by the objective.

SUCCESS CRITERION.  Spearman(`log q_theta`, RMSD) beats both controls on >= 13/19 cells with
a paired CI excluding zero, ON THE IN-BAND SUBSET a selector actually consumes.

FALSIFIER.  If the mean-field control matches or beats `q_theta`, the variational family
contributes nothing beyond factorised smoothing and section 24 is closed.

ALPHA-AS-TEMPERATURE, DERIVED AND THEN MEASURED.  Sprint 16 established alpha's effect is
reproduced by a classical thermostat with Pearson +0.93/+0.98 -- a CORRELATION over six
points.  Here the claim is tested as an exact statement about the object itself: fit the
temperature `T_eff(alpha)` minimising `KL(q_theta || p ~ exp(-E/T))` and report the residual
KL in bits.  A small residual means `q_theta` IS a Boltzmann law of `E` and alpha is nothing
but its temperature.  Both the plain Boltzmann family and the `p0`-tilted family
`p ~ p0 exp(-E/T)` are fitted, because the VQE starts at `p0` and a fair statement must say
which family it lands in.

RUN:  python -m s17.q_repr run
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from s12 import instrument as I
from s14 import vqe_lib as V
from s14 import vqe_run as R
from s15 import seed as SD
from s16 import qphase_lib as QP
from s17 import q_lib as L
from s17 import q_pareto as P

ALPHAS = (0.02, 0.10, 0.25, 1.00)
SEEDS = (0, 1)
#: the "existing candidate pool" the reweighting arm consumes -- a uniform draw, i.e. the
#: retrieval pool's stand-in on this instrument.  Fixed per (target, seed), shared by every
#: reweighting arm so the comparison is on the weights alone.
POOL = 500
PICK = 75


def _bits(nq):
    a = np.arange(1 << nq, dtype=np.int64)
    return ((a[:, None] >> np.arange(nq - 1, -1, -1)[None, :]) & 1).astype(np.int8)


def _H(p):
    p = np.asarray(p, float)
    m = p > 0
    return float(-(p[m] * np.log2(p[m])).sum())


def boltz(E, T, p0=None):
    z = -np.asarray(E, float) / max(float(T), 1e-12)
    z = z - z.max()
    q = np.exp(z) if p0 is None else np.asarray(p0, float) * np.exp(z)
    s = q.sum()
    return q / s if np.isfinite(s) and s > 0 else np.full(len(z), 1.0 / len(z))


def fit_T_kl(q, E, p0=None, lo=1e-4, hi=1e4, iters=60):
    """`T_eff = argmin_T KL(q || Boltzmann(E, T))`, by golden section on log T.

    KL is in BITS.  This is the exact form of the "alpha is a temperature" claim: if the
    residual is small, `q_theta` is a Boltzmann law of the objective and nothing else.
    """
    q = np.asarray(q, float)
    m = q > 0
    qq = q[m]
    lq = np.log2(qq)

    def kl(T):
        b = boltz(E, T, p0)[m]
        b = np.maximum(b, 1e-300)
        return float((qq * (lq - np.log2(b))).sum())

    a, b = np.log(lo), np.log(hi)
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = kl(np.exp(c)), kl(np.exp(d))
    for _ in range(iters):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = kl(np.exp(c))
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = kl(np.exp(d))
    T = float(np.exp((a + b) / 2.0))
    return T, kl(T)


def match_entropy(E, target_bits, p0=None, iters=80):
    """Solve `T` so the Boltzmann law's entropy equals `target_bits`."""
    def h(T):
        return _H(boltz(E, T, p0))
    lo, hi = 1e-6, 1e6
    if target_bits >= h(hi):
        return hi
    if target_bits <= h(lo):
        return lo
    for _ in range(iters):
        m = np.sqrt(lo * hi)
        if h(m) < target_bits:
            lo = m
        else:
            hi = m
    return float(np.sqrt(lo * hi))


def spear_inband(score, rmsd, E, frac=0.01):
    """Global and IN-BAND Spearman. In-band is the only ranking metric that matters here:
    a selector never sees the bulk, it sees the objective's own best `frac`."""
    order = np.argsort(np.asarray(E, float), kind="mergesort")
    m = max(64, int(frac * len(E)))
    sub = order[:m]
    return (float(V.spearman(np.asarray(score)[sub], np.asarray(rmsd)[sub])),
            float(V.spearman(score, rmsd)))


def run_cell(pdb, seed):
    ins = L.inst(pdb)
    E = ins.hamil()
    rmsd = ins.rmsd                                    # ORACLE
    nq = ins.n_qubits
    bits = _bits(nq)
    p0 = QP.untrained_probs(nq, seed)
    rng = SD.stable_rng(pdb, seed, "repr", salt=L.SALT)
    # TWO candidate pools, because "the retrieval pool" has two defensible analogues on this
    # instrument and the answer must not depend on which one is chosen.
    #   `uniform`   -- uniform over the lattice, i.e. uniform over retrieval combinations
    #   `prior`     -- a draw from the retrieval torsion prior's own Boltzmann law, the
    #                  distribution the pipeline's retrieval stage actually induces
    POOLS = {"uniform": np.unique(rng.integers(0, ins.N, POOL)),
             "prior": np.unique(L.prior_draws(ins, float(ins.prior.std()), POOL, rng))}
    pool = POOLS["uniform"]
    out = {"_meta": {"pdb": pdb, "seed": seed, "n": ins.n, "N": ins.N,
                     "fold": ins.fold, "n_pool": int(pool.size),
                     "rho_E_global": float(V.spearman(E, rmsd)),
                     "H_p0_bits": _H(p0)}}

    # -------- the reweighting arms on the FIXED pool, scored on the (M, D) plane ------
    # ENTROPY IS MATCHED ON THE POOL, not on the register.  The first draft matched the
    # register entropy and the realised pool-weight entropies came out 3.30 vs 2.10 bits --
    # a 1.6x diversity confound in the classical control's DISFAVOUR, which is exactly the
    # confound Sprint 16 flagged in `tilt_samples`.  Every control below is scaled to the
    # VQE's own pool-weight entropy by a one-parameter temperature on its log-weight.
    def temper_to(score, H_target, pool):
        s = np.asarray(score, float)[pool]
        s = s - s.max()

        def H_of(b):
            z = b * s
            z = z - z.max()
            w = np.exp(z)
            w /= w.sum()
            return _H(w), w
        if H_of(0.0)[0] <= H_target:            # even uniform is below target: give uniform
            return H_of(0.0)[1]
        lo, hi = 1e-6, 1e6
        if H_of(hi)[0] > H_target:
            return H_of(hi)[1]
        for _ in range(80):
            mid = np.sqrt(lo * hi)
            if H_of(mid)[0] > H_target:
                lo = mid
            else:
                hi = mid
        return H_of(np.sqrt(lo * hi))[1]

    def reweight_point(w, tag, store, pool=pool):
        w = np.asarray(w, float)
        w = w[pool] if w.size == ins.N else w
        w = np.maximum(w, 0.0)
        if w.sum() <= 0:
            w = np.ones_like(w)
        w = w / w.sum()
        acc = {k: [] for k in ("M", "D", "readout", "set_best")}
        for _ in range(6):
            pick = rng.choice(pool, size=PICK, replace=True, p=w)
            r = L.md_plane(ins.ca(pick), ins.nat)
            for k in acc:
                acc[k].append(r[k])
        store[tag] = {k: (float(np.sqrt(np.mean(np.square(v))))
                          if k in ("M", "D", "readout") else float(np.mean(v)))
                      for k, v in acc.items()}
        store[tag]["H_weights_bits"] = _H(w)
        # the argmax-of-weight selection readout, tie-averaged
        mx = w.max()
        tied = pool[w == mx]
        store[tag]["argmax_rmsd"] = float(rmsd[tied].mean())
        store[tag]["argmax_ties"] = int(tied.size)

    store = {}
    reweight_point(np.ones(pool.size) / pool.size, "uniform", store)
    out["reweight_uniform"] = store["uniform"]

    for alpha in ALPHAS:
        key = f"a{alpha}"
        t0 = time.time()
        v = R.run(E, nq, alpha, L.BUDGET, shots=L.SHOTS, ansatz=L.ANSATZ, seed=seed,
                  rmsd=rmsd, bits_per_res=ins.bits_per_res, exact_dist=False,
                  keep_seen=True, lr=0.15, baseline="const")
        an = R.make_ansatz(L.ANSATZ, nq)
        q = np.asarray(an.probs(np.asarray(v["theta"], float)), float)
        q = np.maximum(q, 0.0); q /= q.sum()
        Hq = _H(q)
        logq = np.log(np.maximum(q, 1e-300))

        # ---- ALPHA AS TEMPERATURE: exact KL fits, two families -------------
        T_b, kl_b = fit_T_kl(q, E)
        T_p, kl_p = fit_T_kl(q, E, p0)
        _, kl_p0 = fit_T_kl(q, np.zeros_like(E), p0)   # KL(q || p0) itself, the baseline
        row = {"alpha": alpha, "H_q_bits": Hq,
               "T_eff_boltzmann": T_b, "KL_bits_boltzmann": kl_b,
               "T_eff_p0tilt": T_p, "KL_bits_p0tilt": kl_p,
               "KL_bits_to_p0": kl_p0,
               "mode_rmsd": float(rmsd[int(np.argmax(q))]),
               "mean_rmsd_under_q": float(q @ rmsd),
               "E_pct_under_q": float(q @ QP.ordinal_pct(E)),
               "wall_s": float(time.time() - t0)}

        # ---- P4: the three rankers, global and IN-BAND ----------------------
        Tm = match_entropy(E, Hq)
        qmf, lpmf = L.meanfield_boltzmann(ins, E, Tm)
        Hmf = _H(np.exp(lpmf - lpmf.max()) / np.exp(lpmf - lpmf.max()).sum())
        for tag, score in (("negE", -E), ("logq_theta", logq), ("logq_meanfield", lpmf)):
            ib, gl = spear_inband(score, rmsd, E)
            row[f"rho_inband_{tag}"] = ib
            row[f"rho_global_{tag}"] = gl
        row["H_meanfield_bits"] = Hmf
        row["T_meanfield"] = Tm
        # does q_theta differ from a monotone function of E at all?
        row["rho_logq_vs_negE"] = float(V.spearman(logq, -E))
        row["rho_logqmf_vs_negE"] = float(V.spearman(lpmf, -E))

        # ---- section 24: the reweighting readout on the SAME pool -----------
        # Every control is tempered to the VQE's OWN pool-weight entropy, so the arms differ
        # in WHERE they put the weight and not in HOW MUCH they concentrate it.
        st = {}
        logp0 = np.log(np.maximum(p0, 1e-300))
        for pname, pl in POOLS.items():
            sfx = "" if pname == "uniform" else "_priorpool"
            wq = q[pl] / q[pl].sum() if q[pl].sum() > 0 else np.ones(pl.size) / pl.size
            Ht = _H(wq)
            reweight_point(wq, "vqe" + sfx, st, pool=pl)
            reweight_point(temper_to(-E, Ht, pl), "boltzmann" + sfx, st, pool=pl)
            reweight_point(temper_to(lpmf, Ht, pl), "meanfield" + sfx, st, pool=pl)
            reweight_point(temper_to(logp0 - E / max(T_p, 1e-9), Ht, pl),
                           "p0tilt" + sfx, st, pool=pl)
            # THE MANDATORY MATCHED-RANDOM CONTROL: the VQE's own weight vector, PERMUTED
            # across the pool.  Identical entropy, identical weight multiset, zero
            # information about which candidate is which.
            reweight_point(wq[rng.permutation(pl.size)], "random_matched" + sfx, st,
                           pool=pl)
            reweight_point(np.ones(pl.size) / pl.size, "uniform_control" + sfx, st, pool=pl)
            st["vqe" + sfx]["H_target_bits"] = float(Ht)
        row["reweight"] = st
        out[key] = row
    return out


def run(targets=L.TARGETS9, seeds=SEEDS):
    done = L.ck_load("repr")
    for pdb in targets:
        for sd in seeds:
            key = f"{pdb}_{sd}"
            if key in done:
                print(f"  skip {key}", flush=True)
                continue
            L.gate(key)
            t0 = time.time()
            r = run_cell(pdb, sd)
            L.ck("repr", key, r)
            done[key] = r
            a = r["a0.25"]
            print(f"  {key} {time.time()-t0:.0f}s  KL(q||Boltz) {a['KL_bits_boltzmann']:.3f} "
                  f"bits  KL(q||p0tilt) {a['KL_bits_p0tilt']:.3f}  "
                  f"rho_ib q {a['rho_inband_logq_theta']:+.3f} mf "
                  f"{a['rho_inband_logq_meanfield']:+.3f} E {a['rho_inband_negE']:+.3f}",
                  flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "run":
        tg = sys.argv[2].split(",") if len(sys.argv) > 2 else L.TARGETS9
        run(tg)
    else:
        raise SystemExit(mode)
