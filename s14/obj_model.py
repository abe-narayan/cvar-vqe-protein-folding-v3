"""SPRINT 14, OBJ -- the learned structural objective: features, model, callable.

INFORMATION CONTRACT (this is the whole point of the module).
At inference the scorer sees ONLY:
    * the target sequence,
    * the integer torsion-state configuration (B, n),
    * the target's k-state torsion library (built with the target and everything above
      0.6 identity to it REMOVED -- `torsion_lib2.library_for(seq, k, seq)`),
    * the leakage-safe 1-local empirical prior over those states,
    * model weights trained on OTHER folds.
It never touches native coordinates, native torsions, native distances or RMSD.
RMSD enters only as a TRAINING LABEL and as a post-hoc evaluation axis.

THE MODEL FAMILY.  A sequence-conditioned, distance-binned, many-body pair potential

    E(S) = (1/n) * [ sum_{j-i>=2} w_pair[sep(i,j), cls(i), cls(j), dbin(d_ij)]
                   + sum_i       w_loc[cls(i), rama_cell(phi_i, psi_i)]
                   + w_glob . g(S) ]

It is linear in indicator counts, so ridge fits it exactly from streamed normal
equations, it is SE(3)-invariant by construction (CA distances only), and its Pauli
structure is exactly readable off the sprint-13 locality theorem: a pair term at
sequence separation s = j - i is supported on exactly the s-1 residues strictly between
i and j, i.e. 2*(s-1) qubits at k=4.  Capping the separation caps the Pauli weight.

    from s14.obj_model import LearnedObjective
    obj = LearnedObjective.load("s14/results/obj_lfo_fold3.npz")
    e = obj.score("1CS9", S)          # S: (B, n) int -> (B,) float, lower = better
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I               # noqa: E402
from s13 import qarch_lib as Q                # noqa: E402

# ---------------------------------------------------------------- feature schema
AA_CLASS = {}
for _c, _aas in enumerate(["AVLIMC", "FWY", "STNQ", "KRH", "DE", "G", "P"]):
    for _a in _aas:
        AA_CLASS[_a] = _c
NCLS = 7
CLSPAIR = np.full((NCLS, NCLS), -1, int)
_p = 0
for _i in range(NCLS):
    for _j in range(_i, NCLS):
        CLSPAIR[_i, _j] = CLSPAIR[_j, _i] = _p
        _p += 1
NCLSPAIR = _p                                   # 28

DEDGES = np.array([4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0,
                   9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 20.0])
NDB = len(DEDGES) + 1                           # 17

SEPCLS = np.array([0, 0, 0, 1, 2, 3, 4, 5, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7])
NSEP = 8                                        # sep 2,3,4,5,6,7,{8,9},{10+}

NRAMA = 36                                      # 6 x 6 (phi, psi) cells
NPAIRF = NSEP * NCLSPAIR * NDB                  # 3808
NLOCF = NCLS * NRAMA                            # 252
NGLOB = 8
# composition x global interactions: the ONLY route by which a linear model can make the
# SIGN of its compactness response depend on the sequence, which section 2.3 shows is the
# single thing that decides whether it helps or hurts a target.
NINT = NCLS * NGLOB                             # 56
DIM = NPAIRF + NLOCF + NGLOB + NINT + 1         # + intercept

GLOB_NAMES = ["rg", "e2e_per_n", "contacts8_per_n", "close45_per_n",
              "prior_per_res", "frac_alpha", "frac_beta", "mean_logd"]


def aa_classes(seq):
    return np.array([AA_CLASS.get(c, 0) for c in seq], int)


def rama_cell(phi, psi):
    """(B, n) radians -> (B, n) cell index in [0, 36)."""
    a = np.clip(((phi + np.pi) / (2 * np.pi) * 6).astype(int), 0, 5)
    b = np.clip(((psi + np.pi) / (2 * np.pi) * 6).astype(int), 0, 5)
    return a * 6 + b


class TargetFeatures:
    """Per-target constants: the pair list and its class/separation feature bases."""

    def __init__(self, seq, n=None, max_sep=None):
        self.seq = seq
        self.n = int(n or len(seq))
        cls = aa_classes(seq)[:self.n]
        self.cls = cls
        ii, jj = np.triu_indices(self.n, k=2)
        sep = jj - ii
        if max_sep is not None:
            m = sep <= max_sep
            ii, jj, sep = ii[m], jj[m], sep[m]
        self.ii, self.jj, self.sep = ii, jj, sep
        sc = SEPCLS[np.clip(sep, 0, len(SEPCLS) - 1)]
        cp = CLSPAIR[cls[ii], cls[jj]]
        self.pair_base = ((sc * NCLSPAIR) + cp) * NDB           # (P,)
        self.loc_base = cls * NRAMA + NPAIRF                    # (n,)
        self.npair = len(ii)


class Featurizer:
    """Turns (target, state configuration) into the model's feature representation.

    Holds one `s13.qarch_lib.Space` per target (the held-out torsion library) and its
    leakage-safe empirical prior.  Both are inference-legal.
    """

    def __init__(self, k=4, max_sep=None, cache_dir=None, mask_terminal=True):
        """`mask_terminal` drops the LOCAL (Ramachandran) features of residues 0 and
        n-1 and excludes them from the per-residue prior scalar.

        Justification, verified bit-exactly on the enumerations here (rtol=0, atol=0, on
        7 targets, with residues 1 and n-2 as live controls): under the ideal-geometry
        builder, CA-RMSD is COMPLETELY INVARIANT to the torsion states of residues 0 and
        n-1. Four torsions are inert, not three, so `2*log2(k)` qubits per chain are dead
        (14 of 18 live at n=9, k=4). Those two residues' local features are therefore pure
        noise with respect to the label and only cost capacity.

        Pair features are NOT masked: d(0, j) is a genuine geometric observable that
        depends on residues 1..j-1, so pairs anchored at a terminus are informative even
        though the terminus's own state is not.
        """
        self.k = k
        self.max_sep = max_sep
        self._sp = {}
        self._tf = {}
        self._prior = {}
        self.cache_dir = cache_dir
        self.mask_terminal = bool(mask_terminal)

    def space(self, pdb_id):
        if pdb_id not in self._sp:
            self._sp[pdb_id] = Q.Space(pdb_id, self.k)
        return self._sp[pdb_id]

    def tf(self, pdb_id):
        if pdb_id not in self._tf:
            sp = self.space(pdb_id)
            self._tf[pdb_id] = TargetFeatures(sp.seq, sp.n, self.max_sep)
        return self._tf[pdb_id]

    def prior_table(self, pdb_id):
        """(n, k) leakage-safe occupancy.  Read from the enumeration cache if present."""
        if pdb_id not in self._prior:
            P = None
            if self.cache_dir is not None:
                from s14 import obj_enum as E
                if E.have(pdb_id):
                    z = E.load_enum(pdb_id)
                    if "prior_table" in z.files:
                        P = np.asarray(z["prior_table"], float)
            if P is None:
                P = Q.empirical_prior(self.space(pdb_id))
            self._prior[pdb_id] = P
        return self._prior[pdb_id]

    # -------------------------------------------------------------- geometry
    def geometry(self, pdb_id, S):
        sp = self.space(pdb_id)
        phi, psi = sp.angles(S)
        C = I.build_ca(phi, psi)
        return C, phi, psi

    def rows(self, pdb_id, S, sequence_blind=False):
        """(B, n) states -> (indices (B, m) int32, values (B, m) float32) sparse rows.

        Every row has exactly m = npair + n + NGLOB + 1 entries; duplicate indices are
        summed by the CSR constructor, which is what we want for counts.
        """
        S = np.atleast_2d(np.asarray(S, int))
        B = len(S)
        tf = self.tf(pdb_id)
        n = tf.n
        C, phi, psi = self.geometry(pdb_id, S)
        d = np.linalg.norm(C[:, tf.ii, :] - C[:, tf.jj, :], axis=2)       # (B, P)
        db = np.searchsorted(DEDGES, d).astype(np.int32)
        if sequence_blind:
            sc = SEPCLS[np.clip(tf.sep, 0, len(SEPCLS) - 1)]
            base = ((sc * NCLSPAIR) + CLSPAIR[0, 0]) * NDB
            pidx = base[None, :] + db
        else:
            pidx = tf.pair_base[None, :] + db

        rc = rama_cell(phi, psi)                                          # (B, n)
        if sequence_blind:
            lidx = (rc + NPAIRF).astype(np.int32)
        else:
            lidx = (tf.loc_base[None, :] + rc).astype(np.int32)

        # globals
        cen = C - C.mean(1, keepdims=True)
        rg = np.sqrt((cen ** 2).sum(2).mean(1))
        e2e = np.linalg.norm(C[:, 0] - C[:, -1], axis=1)
        far = tf.sep >= 3
        c8 = (d[:, far] < 8.0).sum(1)
        c45 = (d[:, far] < 4.5).sum(1)
        live = np.ones(n, np.float32)
        if self.mask_terminal and n > 2:
            live[0] = live[-1] = 0.0          # inert for CA-RMSD, verified bit-exactly
        Lp = -np.log(self.prior_table(pdb_id))
        pri = (Lp[np.arange(n)[None, :], S] * live[None, :]).sum(1)
        nl = max(live.sum(), 1.0)
        alpha = (((phi > -np.deg2rad(160)) & (phi < 0)
                  & (psi > -np.deg2rad(120)) & (psi < np.deg2rad(30)))
                 * live[None, :]).sum(1) / nl
        beta = (((phi < 0) & ((psi > np.deg2rad(90)) | (psi < -np.deg2rad(150))))
                * live[None, :]).sum(1) / nl
        mld = np.log(d[:, far] + 1e-6).mean(1) if far.any() else np.zeros(B)
        G = np.stack([rg / 10.0, e2e / (n * 3.8), c8 / n, c45 / n,
                      pri / n, alpha, beta, mld], 1).astype(np.float32)

        gidx = np.arange(NPAIRF + NLOCF, NPAIRF + NLOCF + NGLOB, dtype=np.int32)
        # composition x global interaction block
        frac = np.bincount(tf.cls, minlength=NCLS).astype(np.float32) / n
        if sequence_blind:
            frac = np.zeros(NCLS, np.float32)
        GI = (frac[:, None] * G[:, None, :]).reshape(B, NINT)   # (B, NCLS*NGLOB)
        iidx = np.arange(NPAIRF + NLOCF + NGLOB,
                         NPAIRF + NLOCF + NGLOB + NINT, dtype=np.int32)
        idx = np.concatenate([pidx.astype(np.int32), lidx,
                              np.broadcast_to(gidx, (B, NGLOB)),
                              np.broadcast_to(iidx, (B, NINT)),
                              np.full((B, 1), DIM - 1, np.int32)], 1)
        inv = np.float32(1.0 / n)
        val = np.concatenate([np.full((B, tf.npair), inv, np.float32),
                              np.broadcast_to(live * inv, (B, n)),
                              G, GI, np.ones((B, 1), np.float32)], 1)
        return idx, val

    def csr(self, pdb_id, S, sequence_blind=False):
        from scipy import sparse
        idx, val = self.rows(pdb_id, S, sequence_blind)
        B, m = idx.shape
        indptr = np.arange(0, (B + 1) * m, m, dtype=np.int64)
        X = sparse.csr_matrix((val.ravel(), idx.ravel(), indptr), shape=(B, DIM))
        X.sum_duplicates()
        return X

    def score_with(self, w, pdb_id, S, sequence_blind=False, chunk=100_000):
        """Evaluate the linear objective without materialising a design matrix."""
        S = np.atleast_2d(np.asarray(S, int))
        out = np.empty(len(S), np.float64)
        w = np.asarray(w, np.float64)
        for a in range(0, len(S), chunk):
            idx, val = self.rows(pdb_id, S[a:a + chunk], sequence_blind)
            out[a:a + chunk] = (w[idx] * val).sum(1)
        return out


# ---------------------------------------------------------------- ridge fitting
class RidgeAccumulator:
    """Streaming normal equations for the linear objective."""

    def __init__(self, dim=DIM):
        self.XtX = np.zeros((dim, dim), np.float64)
        self.Xty = np.zeros(dim, np.float64)
        self.n = 0

    def add(self, X, y, center=False, weight=None):
        """Accumulate one target's block.

        `center=True` subtracts that target's own column means and label mean before
        accumulating, EXACTLY and without densifying:

            Xc = X - 1 mu^T,  yc = y - ybar
            Xc^T Xc = X^T X - N mu mu^T
            Xc^T yc = X^T y - N mu ybar

        This turns the fit from "predict absolute structural quality" (which forces the
        model to explain between-target offsets it cannot see) into "predict WITHIN-
        target deviation", which is the ranking objective the search actually consumes.
        """
        y = np.asarray(y, np.float64)
        if weight is not None:
            wsr = np.sqrt(np.asarray(weight, np.float64))
            from scipy import sparse
            X = sparse.diags(wsr) @ X
            y = y * wsr
        N = X.shape[0]
        XtX = (X.T @ X).toarray()
        Xty = X.T @ y
        if center:
            mu = np.asarray(X.sum(0)).ravel() / N
            ybar = float(y.mean())
            XtX = XtX - N * np.outer(mu, mu)
            Xty = Xty - N * mu * ybar
        self.XtX += XtX
        self.Xty += Xty
        self.n += N

    def solve(self, lam=1.0, standardize=True, eps=1e-12):
        """Ridge.  `standardize=True` puts every column on unit within-target variance
        first, so one lambda penalises all feature families comparably.

        This matters here: the count columns and the global scalars differ in variance by
        orders of magnitude, and un-standardised ridge lets the highest-variance columns
        dominate ANY weight vector -- including one fitted to permuted labels, which is
        how a label-permuted null can score like the real model.
        """
        N = max(self.n, 1)
        A = self.XtX / N
        b = self.Xty / N
        if standardize:
            s = np.sqrt(np.maximum(np.diag(A), 0.0))
            s[s < 1e-9] = 1.0
            D = 1.0 / s
            A = A * D[:, None] * D[None, :]
            b = b * D
        A = A + lam * np.eye(len(A))
        if not standardize:
            A[-1, -1] -= lam                    # do not penalise the intercept
        # With per-target centering the intercept column is exactly zero, so exempting
        # it would leave a zero row AND column and make the system singular. It only
        # shifts every score by a constant, so rankings are unaffected either way.
        w = np.linalg.solve(A, b)
        return w * D if standardize else w

    def solve_path(self, lams, standardize=True):
        return {float(l): self.solve(l, standardize) for l in lams}


class LearnedObjective:
    """The deliverable callable.  Deterministic, importable, inference-legal."""

    def __init__(self, w, k=4, max_sep=None, sequence_blind=False, meta=None,
                 mask_terminal=True):
        self.w = np.asarray(w, np.float64)
        self.k = int(k)
        self.max_sep = max_sep
        self.sequence_blind = bool(sequence_blind)
        self.mask_terminal = bool(mask_terminal)
        self.meta = dict(meta or {})
        self.fz = Featurizer(k=self.k, max_sep=self.max_sep, cache_dir=True,
                             mask_terminal=self.mask_terminal)

    def score(self, pdb_id, S):
        """(B, n) integer torsion states -> (B,) energy-like scores, lower is better."""
        return self.fz.score_with(self.w, pdb_id, S, self.sequence_blind)

    __call__ = score

    def save(self, path):
        np.savez(path, w=self.w, k=self.k,
                 max_sep=-1 if self.max_sep is None else self.max_sep,
                 sequence_blind=self.sequence_blind,
                 mask_terminal=self.mask_terminal,
                 meta=np.array([repr(self.meta)], dtype=object), allow_pickle=True)
        return path

    @classmethod
    def load(cls, path):
        z = np.load(path, allow_pickle=True)
        ms = int(z["max_sep"])
        return cls(z["w"], int(z["k"]), None if ms < 0 else ms,
                   bool(z["sequence_blind"]),
                   mask_terminal=bool(z["mask_terminal"]))

    # ---- Pauli locality, from the sprint-13 exact locality theorem -----------
    def pauli_support(self, n):
        """Qubit support of each retained pair term at k=4 (2 qubits per residue).

        d_ij depends on exactly the j-i-1 residues strictly between i and j, so a
        separation-s pair term touches 2*(s-1) qubits.  Returns {sep: n_qubits}.
        """
        ms = self.max_sep or (n - 1)
        return {int(s): int(2 * (s - 1)) for s in range(2, min(ms, n - 1) + 1)}
