"""SPRINT 19 / AGENT B -- machinery for the CONTINUOUS-TORSION CVaR-VQE sampler test.

Pre-registered in `s19/PREREG_B.md` before this file existed.

WHAT IS NEW HERE, AND WHY IT IS NOT A REPEAT OF SPRINT 17/18
-----------------------------------------------------------
Every previous quantum arm in this programme optimised or sampled on the **k=4 lattice
register**: a discrete space of `k^n` configurations, indexed by a bitstring, with an
ARBITRARY 2-bit code naming each torsion state.  Sprint 18's headline was that the one
interval in the degree-1 story which excluded zero measured *the arbitrariness of that
encoding*, at gauge percentile 0.00.

This module does not use that register at all.  The configuration is

    z = (phi_1..phi_n, psi_1..psi_n)  in  T^{2n}          -- CONTINUOUS

exactly the object `s15/align_lib.fit` optimises, and the objective is the deployed
distogram functional evaluated on it.  The circuit's bitstring is a **latent that selects
which conformer basin each residue occupies**; the angle itself is then drawn CONTINUOUSLY
from a von Mises component with finite concentration, so the sampler's law is a continuous
density with FULL SUPPORT on the torus:

    p_theta(phi, psi) = sum_b p_theta(b) prod_i q_{i,b_i}(phi_i, psi_i)

The only discrete object is `b`, and its only arbitrariness is `b_i <-> 1-b_i` per residue --
a Z2^n gauge, which `qb_gauge` tests directly.  That is a 2-element per-residue orbit against
Sprint 18's 24, and it is the *whole* labelling freedom, not a sample of it.

WHY A BASIN LATENT.  What an entangled latent can supply that a product latent cannot is
CORRELATION BETWEEN WHICH BASIN NEIGHBOURING RESIDUES OCCUPY.  That is precisely the
structure `s19/BRIEF.md` section 3 names: "a whole region flipping between two conformer
families moves many pairs the same way at once."  If a structured quantum sampler is worth
anything on this problem, this is where it shows.  It is built to make its failure mean
something.

BUDGET ACCOUNTING (Sprint 14's rule, unchanged).  One objective evaluation = one read of
`E` at one continuous configuration.  Circuit simulation, latent sampling, grad-log-p and
the von Mises draws are FREE.  Every budgeted arm gets B = 8192.

ORACLE LABELLING.  Every CA-RMSD, every `M`, every coverage count and every generation
ceiling reads the native and is post-hoc scoring only.  `D`, `n_distinct`, entropy, the
objective, the basin fits and every selection decision are NATIVE-FREE.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(HERE, "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

from s12 import instrument as I                 # noqa: E402
from s15 import distcal as C                    # noqa: E402
from s15 import seed as SD                      # noqa: E402
from s14 import retprior as RP                  # noqa: E402

SALT = "s19qb"
BUDGET = 8192
SHOTS = 512
LR = 0.15
M_SEL = 75


# ============================================================ checkpointing
_CK = {}


def ck(tag, key, value):
    """Merge-on-write incremental checkpoint into `s19/results/qb_<tag>.json`."""
    path = os.path.join(RESULTS, f"qb_{tag}.json")
    if os.path.exists(path):
        try:
            with open(path) as fh:
                disk = json.load(fh)
            mem = _CK.setdefault(tag, {})
            for k, v in disk.items():
                mem.setdefault(k, v)
        except Exception:
            pass
    d = _CK.setdefault(tag, {})
    d[key] = value
    d["_written"] = time.strftime("%Y-%m-%d %H:%M:%S")
    tmp = path + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(d, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    for _ in range(30):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            time.sleep(0.4)
    with open(path, "w") as fh:
        json.dump(d, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))


def load(tag):
    path = os.path.join(RESULTS, f"qb_{tag}.json")
    if os.path.exists(path):
        with open(path) as fh:
            d = json.load(fh)
        _CK[tag] = d
        return d
    return _CK.setdefault(tag, {})


# ============================================================ the instrument
_DATA = None
_DEB = None


def gather_full():
    """`C.gather` + the leave-fold-out separation debias fitted on the FULL 126 (subset trap)."""
    global _DATA, _DEB
    if _DATA is None:
        tg = I.targets()
        pdbs = [t["pdb"] for t in tg]
        _DATA = C.gather(tg)
        deb = {}
        for f in sorted({int(t["fold"]) for t in tg}):
            train = [p for p in pdbs if _DATA[p]["fold"] != f]
            fn, _ = C.fit_correction(_DATA, train, "sep")
            deb[f] = (lambda sp, _fn=fn: -_fn(np.zeros_like(np.asarray(sp, float)), sp))
        _DEB = deb
    return _DATA, _DEB


def wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


# ------------------------------------------------------ 2-component von Mises basins
def fit_basins(PHI, PSI, rng, kappa_cap=40.0, kappa_floor=1.0):
    """NATIVE-FREE.  Per residue, a 2-component von Mises mixture over (phi, psi).

    `PHI`,`PSI` are the (m, n) retrieval-pool torsions.  Circular k-means with k=2 on the
    4-D embedding (cos/sin of both angles), then moment-matched concentrations per
    component per angle.  `kappa` is capped ABOVE (so no component becomes a delta and the
    mixture keeps full support on the torus) and BELOW (so a component is not uniform).

    Returns mu (n, 2, 2), kap (n, 2, 2), w (n, 2) -- the marginal basin occupancy.
    """
    PHI = wrap(PHI); PSI = wrap(PSI)
    m, n = PHI.shape
    mu = np.zeros((n, 2, 2))
    kap = np.zeros((n, 2, 2))
    w = np.zeros((n, 2))
    for i in range(n):
        X = np.column_stack([np.cos(PHI[:, i]), np.sin(PHI[:, i]),
                             np.cos(PSI[:, i]), np.sin(PSI[:, i])])
        # circular k-means, k=2, deterministic init: the two most distant rows
        d0 = ((X - X.mean(0)) ** 2).sum(1)
        a = int(np.argmax(d0))
        b = int(np.argmax(((X - X[a]) ** 2).sum(1)))
        cen = X[[a, b]].copy()
        lab = np.zeros(m, int)
        for _ in range(40):
            dd = ((X[:, None, :] - cen[None]) ** 2).sum(2)
            nl = np.argmin(dd, 1)
            if (nl == lab).all():
                break
            lab = nl
            for c in (0, 1):
                if (lab == c).any():
                    cen[c] = X[lab == c].mean(0)
                else:                                  # empty cluster: reseed on the far row
                    cen[c] = X[int(np.argmax(((X - cen[1 - c]) ** 2).sum(1)))]
        for c in (0, 1):
            sel = lab == c
            if sel.sum() < 2:
                sel = np.ones(m, bool)
            w[i, c] = float(max(sel.sum(), 1)) / m
            for t, A in ((0, PHI[:, i]), (1, PSI[:, i])):
                S, Cc = np.sin(A[sel]).mean(), np.cos(A[sel]).mean()
                mu[i, c, t] = np.arctan2(S, Cc)
                R = float(np.hypot(S, Cc))
                R = min(max(R, 1e-6), 0.999)
                # standard von Mises A^{-1}(R) approximation (Banerjee et al.)
                k = R * (2.0 - R * R) / (1.0 - R * R)
                kap[i, c, t] = float(np.clip(k, kappa_floor, kappa_cap))
        w[i] /= w[i].sum()
    return mu, kap, w


def draw_from_basins(bits, mu, kap, rng):
    """bits (B, n) uint8 -> (phi, psi) each (B, n), CONTINUOUS von Mises draws."""
    bits = np.asarray(bits, np.int64)
    B, n = bits.shape
    ii = np.arange(n)[None, :]
    phi = rng.vonmises(mu[ii, bits, 0], kap[ii, bits, 0])
    psi = rng.vonmises(mu[ii, bits, 1], kap[ii, bits, 1])
    return wrap(phi), wrap(psi)


# ============================================================ the objective
class Obj:
    """The DEPLOYED distogram objective on CONTINUOUS torsions, with a hard eval budget.

    `E(phi, psi) = sum_p (d_p - dhat_p)^2 / sd_p^2` -- identical to `s15/align_lib.fit`'s
    `f` at kind="squared", wpair=None.  Every call increments `used` by the batch size and
    truncates at `budget`; `seen_phi/seen_psi` keep the emitted set for the readout.
    """

    def __init__(self, dhat, sd, i, j, budget=BUDGET, keep=True):
        self.dhat = np.asarray(dhat, float)
        self.inv2 = 1.0 / np.asarray(sd, float) ** 2
        self.i = np.asarray(i, int)
        self.j = np.asarray(j, int)
        self.budget = int(budget)
        self.used = 0
        self.keep = keep
        self._phi, self._psi, self._e = [], [], []

    @property
    def left(self):
        return max(0, self.budget - self.used)

    def raw(self, phi, psi):
        """Unbudgeted, unrecorded evaluation.  For post-hoc scoring only."""
        from core import project as pj
        phi = np.atleast_2d(np.asarray(phi, float))
        psi = np.atleast_2d(np.asarray(psi, float))
        CA = np.asarray(pj.build_ca_exact(phi, psi), float)
        rv = CA[:, self.i, :] - CA[:, self.j, :]
        d = np.sqrt((rv * rv).sum(-1))
        return ((d - self.dhat) ** 2 * self.inv2).sum(-1), CA

    def __call__(self, phi, psi):
        """Budgeted.  Returns the energies of as many rows as the budget allowed."""
        phi = np.atleast_2d(np.asarray(phi, float))
        psi = np.atleast_2d(np.asarray(psi, float))
        k = min(len(phi), self.left)
        if k <= 0:
            return np.zeros(0)
        phi, psi = phi[:k], psi[:k]
        e, _CA = self.raw(phi, psi)
        self.used += k
        if self.keep:
            self._phi.append(phi.copy()); self._psi.append(psi.copy()); self._e.append(e)
        return e

    def seen(self):
        if not self._phi:
            return np.zeros((0, 0)), np.zeros((0, 0)), np.zeros(0)
        return (np.concatenate(self._phi), np.concatenate(self._psi),
                np.concatenate(self._e))


# ============================================================ the readout
def readout(phi, psi, e, tgt, m=M_SEL, tag=""):
    """THE FOUR OUTCOMES, SEPARATELY (brief section 6).  ORACLE columns are post-hoc.

    generation ceiling (set best) | (M, D) plane | selection ceiling | realised.
    Selection is the FROZEN NATIVE-FREE shipped Bayes-risk selector on the emitted set.
    Realised applies the frozen terminal operator: coordinate average of the selected m,
    then the ideal-geometry projection (the repair operator).
    """
    from s17 import q_lib as QL
    from core import project as pj
    nat = tgt["nat"]
    phi = np.asarray(phi, float); psi = np.asarray(psi, float)
    if len(phi) == 0:
        return {"tag": tag, "n_emitted": 0}
    CA = np.asarray(pj.build_ca_exact(phi, psi), float)
    rr = I.kabsch_rmsd_batch(CA, nat)                                   # ORACLE
    # --- native-free selection: the shipped Bayes-risk score on the emitted set
    D = I.pair_dists(CA, tgt["i"], tgt["j"])
    sc = I.shipped_score(tgt["dg"], D)
    mm = min(m, len(phi))
    sel = np.argsort(sc, kind="stable")[:mm]
    # THE MIN-OF-N NULL'S BAND, NAMED.  A sampler emitting 8192 configurations has a
    # min-of-N advantage over the 500-window incumbent pool that has nothing to do with
    # the sampler.  `sub500` is the identical statistic on a stable 500-subsample, so the
    # incumbent comparison can be made at matched count.
    sub = SD.stable_rng(tgt["pdb"], tag, "sub500", salt=SALT).permutation(len(phi))[:500]
    Wsel = CA[sel]
    md = QL.md_plane(Wsel, nat)
    avg, _b = I.coordinate_average(Wsel)
    pr = I.project(np.asarray(avg, float), tgt["seq"], tgt["fold"])
    out = {
        "tag": tag,
        "n_emitted": int(len(phi)),
        "n_distinct": int(len(phi)),
        "obj_best": float(np.min(e)) if len(e) else float("nan"),
        "obj_mean": float(np.mean(e)) if len(e) else float("nan"),
        # 1. GENERATION -- ORACLE
        "gen_best_ORACLE": float(rr.min()),
        "gen_mean_ORACLE": float(rr.mean()),
        "gen_p05_ORACLE": float(np.percentile(rr, 5)),
        "gen_best_sub500_ORACLE": float(rr[sub].min()),
        "sel_best_sub500_ORACLE": float(
            rr[sub][np.argsort(sc[sub], kind="stable")[:min(m, len(sub))]].min()),
        "cov_2.5_ORACLE": int((rr < 2.5).sum()),
        "cov_2.0_ORACLE": int((rr < 2.0).sum()),
        # 2. the (M, D) plane on the SELECTED set -- M ORACLE, D native-free
        "M": md["M"], "D": md["D"], "readout_MD_ORACLE": md["readout"],
        "identity_resid": md["identity_resid"],
        # 3. SELECTION ceiling -- ORACLE
        "sel_best_ORACLE": float(rr[sel].min()),
        "sel_mean_ORACLE": float(rr[sel].mean()),
        # 4. REALISED -- ORACLE scoring of a native-free decision
        "avg_ORACLE": float(I.ca_rmsd(np.asarray(avg, float), nat)),
        "realised_ORACLE": float(I.ca_rmsd(np.asarray(pr["ca"], float), nat)),
    }
    return out


# ============================================================ target bundle
def target(pdb, budget=BUDGET):
    data, deb = gather_full()
    d = data[pdb]
    fold, n, seq = int(d["fold"]), int(d["n"]), d["seq"]
    dhat = np.maximum(d["dhat"] - deb[fold](d["sep"]), 2.0)
    PHI, PSI, _sim = RP.windows(pdb, "pool")
    rng = SD.stable_rng(pdb, "basins", salt=SALT)
    mu, kap, w = fit_basins(PHI, PSI, rng)
    return {"pdb": pdb, "n": n, "seq": seq, "fold": fold,
            "i": d["i"], "j": d["j"], "sd": d["sd"], "dhat": dhat,
            "nat": np.asarray(d["nat"], float),
            "dg": I.distogram(pdb, seq, fold),
            "PHI": PHI, "PSI": PSI, "mu": mu, "kap": kap, "wmarg": w}


def new_obj(tgt, budget=BUDGET):
    return Obj(tgt["dhat"], tgt["sd"], tgt["i"], tgt["j"], budget=budget)
