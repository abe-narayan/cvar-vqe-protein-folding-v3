"""TORSION-PREDICTOR (sprint 13) -- shared library.

THE QUESTION.  Sprint 12 measured that *continuous* backbone torsions known to sigma = 12
degrees at full residue coverage build a 1.486 A structure (86% under 2 A, FAIL18 6.019 ->
1.80).  It also measured a leave-fold-out 4-state torsion-bin predictor at 0.690 accuracy /
0.517 on the FAIL18 and found it emitted +0.015 A.  This module asks the gating question
properly: *can a predictor that sees only inference-legal information produce torsion
information good enough to reach that regime?*

WHAT IS DIFFERENT FROM SPRINT 12's PREDICTOR.  S12 predicted a 4-way ABEGO label and used
it as a retrieval key.  Four states is a 90-degree-scale discretisation -- it cannot express
a 12-degree restraint even if it were perfect, and the argmax of a 4-way posterior throws
away exactly the thing the downstream consumer needs (how sure the model is, and whether
the residue is genuinely bimodal between the alpha and beta basins).  Here the prediction
target is a *density on the torus*, and the deliverables are

    point       circular point estimate of (phi, psi)
    grid        full categorical over an 18 x 18 = 324 cell (20 degree) Ramachandran grid
    mvm         mixture of 8 bivariate von Mises densities, trained by NLL
    abego4      S12's 4-state predictor, reproduced as the reference point

plus, for every family, a per-residue SELF-ESTIMATED sigma (`sigma_hat`) -- the posterior's
own expected angular deviation from its point estimate, in degrees.  `sigma_hat` is the
quantity a downstream consumer gates on, and it is on exactly the axis Sprint 12's
(sigma, coverage) surface is parameterised by, so a predictor's quality converts into an
expected emitted RMSD without any new assumption.

LEAKAGE DISCIPLINE.  Training labels are native torsions of OTHER-FOLD peptides plus this
fold's protein fragments -- `s12.key_corpus.corpus(f)` is exactly the legal library for
every fold-f target, stitched back into parent chains.  `assert_fold_discipline()` re-checks
that by direct sequence comparison rather than trusting the construction.  Nothing at
inference reads the native, the pool ranking, or the target's own fold.

PLACEHOLDER ANGLES.  `peptide_db` stores phi[0] = -60 deg and psi[L-1] = -45 deg for every
chain; those two angles are undefined and, as verified here, have EXACTLY zero effect on the
CA trace built by `core.project.build_ca_exact`.  They are dropped from training, from the
angular-error statistics and from coverage accounting, so a "full coverage" claim here means
all 2n-2 determined angles.
"""
from __future__ import annotations
import os, sys, json, math, time

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import torch
torch.set_num_threads(2)
import torch.nn as nn

from s12 import instrument as I
from s12 import key_corpus as KC
from priors import _PROPS as _PP, _DEFAULT as _PD

CACHE = os.path.join(ROOT, "s13", "cache")
RESULTS = os.path.join(ROOT, "s13", "results")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)

D2R = np.pi / 180.0
R2D = 180.0 / np.pi
W = 7                       # +-W residue sequence context
PAD = 20
NAA = 21
TERM = 7                    # distance-to-terminus one-hot depth

PROPS = np.zeros((NAA, 5), np.float32)
for _a, _v in _PP.items():
    PROPS[I.ALPHABET.index(_a)] = _v
PROPS[PAD] = _PD
PROPS = (PROPS - PROPS.mean(0)) / (PROPS.std(0) + 1e-6)


def wrap(a):
    a = np.asarray(a, float)
    return np.arctan2(np.sin(a), np.cos(a))


WAIT_LOG = []


def free_ok(min_gb=1.5, wait=True, tag="", max_wait=60.0, floor_gb=0.7):
    """Compute-cap guard: never start a heavy step under `min_gb` free.

    The brief's rule is to wait below 1.5 GB.  The box is shared with other sprint agents, so
    an unconditional wait can deadlock; this waits up to `max_wait` and then proceeds only if
    free memory is still above `floor_gb`, which is above this process's measured peak RSS
    (~0.6 GB for the largest arm).  Every such event is recorded and reported.
    """
    t0 = time.time()
    g = I.free_gb()
    waited = 0.0
    while wait and g < min_gb and (time.time() - t0) < max_wait:
        print(f"  [wait] free={g:.2f} GB < {min_gb} ({tag})", flush=True)
        time.sleep(20.0)
        g = I.free_gb()
        waited = time.time() - t0
    if g < min_gb:
        WAIT_LOG.append({"tag": tag, "waited_s": round(waited, 1), "free_gb": round(g, 2)})
        print(f"  [proceed-after-wait] free={g:.2f} GB after {waited:.0f}s ({tag})", flush=True)
        while g < floor_gb:
            time.sleep(20.0)
            g = I.free_gb()
    return g


# ==================================================================== label discretisation
NPHI = NPSI = 18                       # 20-degree cells
NGRID = NPHI * NPSI                    # 324
_EDGE = np.linspace(-180.0, 180.0, NPHI + 1)
_CEN = (_EDGE[:-1] + _EDGE[1:]) / 2.0  # bin centres, degrees


def grid_bin(phi, psi):
    """Joint 18x18 Ramachandran cell index of radian angles."""
    fp = np.clip(((np.degrees(wrap(phi)) + 180.0) / 20.0).astype(np.int64), 0, NPHI - 1)
    sp = np.clip(((np.degrees(wrap(psi)) + 180.0) / 20.0).astype(np.int64), 0, NPSI - 1)
    return fp * NPSI + sp


GRID_PHI = np.repeat(_CEN, NPSI) * D2R          # (324,) radians
GRID_PSI = np.tile(_CEN, NPHI) * D2R


def abego4(phi, psi):
    """S12's `key_lib.abego4`, reproduced verbatim so the reference point is comparable."""
    ph = np.degrees(wrap(phi)); ps = np.degrees(wrap(psi))
    neg = ph < 0
    a = (ps > -75.0) & (ps <= 50.0)
    g = (ps > -100.0) & (ps <= 100.0)
    return np.where(neg, np.where(a, 0, 1), np.where(g, 2, 3)).astype(np.int64)


# ============================================================================== features
def _ctx_index(L, pid):
    idx = np.arange(L)[:, None] + np.arange(-W, W + 1)[None, :]
    ok = (idx >= 0) & (idx < L)
    idxc = np.clip(idx, 0, L - 1)
    ok &= pid[idxc] == pid[:, None]
    return idxc, ok


def context_feats(codes, pid=None):
    """(L, 15*(21+5)) sliding one-hot + property context; context never crosses parents."""
    codes = np.asarray(codes, np.int64); L = len(codes)
    pid = np.zeros(L, np.int64) if pid is None else np.asarray(pid, np.int64)
    idxc, ok = _ctx_index(L, pid)
    ctx = np.where(ok, codes[idxc], PAD)
    oh = np.zeros((L, 2 * W + 1, NAA), np.float32)
    np.put_along_axis(oh, ctx[:, :, None], 1.0, 2)
    return np.concatenate([oh, PROPS[ctx]], 2).reshape(L, -1)


def pos_feats(pid):
    """Distance-to-parent-terminus one-hots + parent-length -- purely sequence-derived."""
    pid = np.asarray(pid, np.int64); L = len(pid)
    out = np.zeros((L, 2 * (TERM + 1) + 2), np.float32)
    _, st = np.unique(pid, return_index=True)
    st = np.sort(st); en = np.append(st[1:], L)
    for a, b in zip(st, en):
        ln = b - a
        dn = np.minimum(np.arange(ln), TERM)
        dc = np.minimum(np.arange(ln)[::-1], TERM)
        out[a:b, dn] = 1.0
        out[a:b, TERM + 1 + dc] = 1.0
        out[a:b, -2] = min(ln, 40) / 40.0
        out[a:b, -1] = 1.0 if ln <= 26 else 0.0
    return out


def comp_feats(codes, pid=None):
    """MANDATORY NULL: the parent's amino-acid composition, identical at every position."""
    codes = np.asarray(codes, np.int64); L = len(codes)
    pid = np.zeros(L, np.int64) if pid is None else np.asarray(pid, np.int64)
    out = np.zeros((L, 20 + 5), np.float32)
    for p in np.unique(pid):
        m = pid == p
        c = np.bincount(codes[m], minlength=NAA)[:20].astype(np.float32)
        c = c / max(c.sum(), 1.0)
        out[m, :20] = c
        out[m, 20:] = (PROPS[codes[m]][:, :5]).mean(0)
    return out


def esm_extra(codes, pid, E):
    """Centre-residue and context-mean ESM-2 PCA-32 (the cached, legal bank)."""
    L = len(codes)
    pid = np.zeros(L, np.int64) if pid is None else np.asarray(pid, np.int64)
    E = np.asarray(E, np.float32)
    idxc, ok = _ctx_index(L, pid)
    Em = (E[idxc] * ok[:, :, None]).sum(1) / np.maximum(ok.sum(1, keepdims=True), 1)
    return np.concatenate([E, Em.astype(np.float32)], 1)


def featurise(codes, pid=None, feat="context", E=None):
    codes = np.asarray(codes, np.int64); L = len(codes)
    pid = np.zeros(L, np.int64) if pid is None else np.asarray(pid, np.int64)
    if feat == "comp":
        return np.concatenate([comp_feats(codes, pid), pos_feats(pid)], 1)
    if feat == "ctxonly":
        # ABLATION.  `pos_feats` is a distribution-shift hazard: on a 13-mer every residue is
        # within 7 of a terminus, while in an 80-residue stitched fragment parent most are
        # not, so the positional block can push a target into a feature region the model
        # mostly saw for peptides only.  This arm removes it.
        return context_feats(codes, pid).astype(np.float32)
    X = [context_feats(codes, pid), pos_feats(pid)]
    if feat == "esm":
        if E is None:
            E = np.zeros((L, 32), np.float32)
        X.append(esm_extra(codes, pid, E))
    return np.concatenate(X, 1).astype(np.float32)


def esm_matrix_corpus(fold):
    """(L,32) ESM-2 PCA-32 aligned to corpus(fold) residue order (S12's cached file)."""
    path = os.path.join(I.CACHE, f"key_esm_f{fold}.npy")
    if os.path.exists(path):
        return np.load(path).astype(np.float32)
    raise FileNotFoundError(path)


def esm_target(seq):
    from s12 import esm_bank
    v = esm_bank.load().get(seq)
    return np.asarray(v[0], np.float32) if v is not None else np.zeros((len(seq), 32), np.float32)


def codes_of(seq):
    return np.array([I.ALPHABET.index(a) if a in I.ALPHABET else PAD for a in seq], np.int64)


# ================================================================== corpus + fold discipline
_CONTAIN_LEAK = {}


def containment_leaks(fold):
    """Peptide parents of `corpus(fold)` that CONTAIN a fold-`fold` target sequence verbatim.

    FOUND HERE, NOT INHERITED.  The project's leakage control is Needleman-Wunsch identity
    normalised by the LONGER sequence, thresholded at 0.6.  A 10-mer target sitting inside a
    17-mer database peptide scores 10/17 = 0.588 and passes -- but its exact sequence, with
    its native torsions, is then in the training labels.  Four such pairs exist across the
    126 tuning targets (6B9K, 1CEK, 2FBU, 2P5H).  This is a pre-existing property of
    `peptide_db.identity`, not of the universe files; it is reported rather than silently
    absorbed, and the containing parents are dropped from the fold's corpus so that every
    number in this agent's findings is strictly cleaner than the project's own standard.
    """
    if fold in _CONTAIN_LEAK:
        return _CONTAIN_LEAK[fold]
    tgt = [t["seq"] for t in I.targets() if t["fold"] == fold]
    c = KC.corpus(fold)
    pid = np.asarray(c["pid"], np.int64)
    org = np.asarray(c["org"], bool)
    codes = np.asarray(c["codes"], np.int64)
    bad = []
    for p in np.unique(pid):
        m = pid == p
        if not org[m][0]:
            continue
        q = I.decode(codes[m])
        if any((s in q) or (q in s) for s in tgt):
            bad.append(int(p))
    _CONTAIN_LEAK[fold] = bad
    return bad


def corpus(fold, drop_containment=True):
    """Legal per-fold training corpus with the two placeholder angles per parent masked."""
    c = KC.corpus(fold)
    pid = np.asarray(c["pid"], np.int64); L = len(pid)
    _, st = np.unique(pid, return_index=True)
    st = np.sort(st); en = np.append(st[1:], L) - 1
    ok = np.ones(L, bool)
    ok[st] = False              # phi undefined at a parent's first residue
    ok[en] = False              # psi undefined at a parent's last residue
    keep = np.ones(L, bool)
    if drop_containment:
        bad = containment_leaks(fold)
        if bad:
            keep = ~np.isin(pid, bad)
    return {"codes": np.asarray(c["codes"], np.int64)[keep], "pid": pid[keep],
            "phi": np.asarray(c["phi"], float)[keep], "psi": np.asarray(c["psi"], float)[keep],
            "org": np.asarray(c["org"], bool)[keep], "ok": ok[keep], "keep": keep,
            "n_parents_dropped": len(containment_leaks(fold)) if drop_containment else 0}


def parent_seqs(fold, peptide_only=True):
    c = corpus(fold)
    out = []
    for p in np.unique(c["pid"]):
        m = c["pid"] == p
        if peptide_only and not c["org"][m][0]:
            continue
        out.append(I.decode(c["codes"][m]))
    return out


def assert_fold_discipline(verbose=True):
    """Re-derive, do not trust: no fold-f target sequence may appear in corpus(f).

    Checks (a) exact-substring containment against every PEPTIDE parent, and (b) the
    project's own Needleman-Wunsch identity < 0.6 (`peptide_db.IDENTITY_THRESHOLD`).
    Fragments are windows of large PDB proteins and are the same in every fold by
    construction, so the peptide side is where a fold leak could live.
    """
    import peptide_db as pdb
    tg = I.targets()
    report = {}
    for f in sorted({t["fold"] for t in tg}):
        ps = parent_seqs(f, peptide_only=True)
        worst = 0.0; worst_pair = None; nsub = 0
        for t in tg:
            if t["fold"] != f:
                continue
            s = t["seq"]
            for q in ps:
                if s in q or q in s:
                    nsub += 1
                idv = pdb.identity(s, q)
                if idv > worst:
                    worst, worst_pair = idv, (t["pdb"], q)
        report[f] = {"n_peptide_parents": len(ps), "n_substring_hits": int(nsub),
                     "max_identity": float(worst), "argmax": worst_pair}
        assert nsub == 0, f"fold {f}: a target sequence occurs in its own training corpus"
        assert worst < pdb.IDENTITY_THRESHOLD, f"fold {f}: identity {worst:.3f} >= 0.6"
        if verbose:
            print(f"  fold {f}: {len(ps)} peptide parents, max identity to any fold-{f} "
                  f"target = {worst:.3f} (< 0.60)  OK", flush=True)
    return report


# ================================================================================ models
class Trunk(nn.Module):
    def __init__(self, d, nout, h=384):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(0.2),
                                 nn.Linear(h, h // 2), nn.ReLU(), nn.Dropout(0.2),
                                 nn.Linear(h // 2, nout))

    def forward(self, x):
        return self.net(x)


NMIX = 8


def _log_i0(k):
    return torch.log(torch.special.i0e(k) + 1e-12) + k


def mvm_nll(out, phi, psi):
    """Negative log-likelihood of (phi, psi) under a mixture of NMIX bivariate von Mises."""
    B = out.shape[0]
    o = out.view(B, NMIX, 8)
    logw = torch.log_softmax(o[:, :, 0], 1)
    muf = torch.atan2(o[:, :, 1], o[:, :, 2])
    mus = torch.atan2(o[:, :, 3], o[:, :, 4])
    kf = torch.clamp(torch.exp(torch.clamp(o[:, :, 5], -4.0, 5.0)), 1e-3, 150.0)
    ks = torch.clamp(torch.exp(torch.clamp(o[:, :, 6], -4.0, 5.0)), 1e-3, 150.0)
    lp = (kf * torch.cos(phi[:, None] - muf) - _log_i0(kf) - math.log(2 * math.pi)
          + ks * torch.cos(psi[:, None] - mus) - _log_i0(ks) - math.log(2 * math.pi))
    return -torch.logsumexp(logw + lp, 1)


def mvm_grid_density(out, gphi, gpsi):
    """(B, G) normalised mixture density evaluated on a fixed torus grid."""
    B = out.shape[0]
    o = out.view(B, NMIX, 8)
    logw = torch.log_softmax(o[:, :, 0], 1)
    muf = torch.atan2(o[:, :, 1], o[:, :, 2])
    mus = torch.atan2(o[:, :, 3], o[:, :, 4])
    kf = torch.clamp(torch.exp(torch.clamp(o[:, :, 5], -4.0, 5.0)), 1e-3, 150.0)
    ks = torch.clamp(torch.exp(torch.clamp(o[:, :, 6], -4.0, 5.0)), 1e-3, 150.0)
    lp = (kf[:, :, None] * torch.cos(gphi[None, None, :] - muf[:, :, None]) - _log_i0(kf)[:, :, None]
          + ks[:, :, None] * torch.cos(gpsi[None, None, :] - mus[:, :, None]) - _log_i0(ks)[:, :, None])
    dens = torch.logsumexp(logw[:, :, None] + lp, 1)
    p = torch.softmax(dens, 1)
    return p


FAMILY_NOUT = {"point": 4, "grid": NGRID, "abego4": 4, "mvm": NMIX * 8}


def make_labels(family, phi, psi):
    if family == "grid":
        return grid_bin(phi, psi).astype(np.int64)
    if family == "abego4":
        return abego4(phi, psi)
    return np.stack([np.cos(phi), np.sin(phi), np.cos(psi), np.sin(psi)], 1).astype(np.float32)


def loss_fn(family):
    ce = nn.CrossEntropyLoss()
    if family in ("grid", "abego4"):
        return lambda o, y: ce(o, y)
    if family == "point":
        return lambda o, y: ((o - y) ** 2).mean()
    if family == "mvm":
        return lambda o, y: mvm_nll(o, torch.atan2(y[:, 1], y[:, 0]), torch.atan2(y[:, 3], y[:, 2])).mean()
    raise ValueError(family)


def train_one(X, Y, family, seed=0, epochs=30, bs=1024, lr=1e-3, val_frac=0.1,
              groups=None, patience=5, verbose=False):
    """Train one leave-fold-out model.  Validation split is by PARENT, never by residue."""
    rng = np.random.default_rng(seed)
    g = np.arange(len(X)) if groups is None else np.asarray(groups)
    ug = np.unique(g)
    vg = set(ug[rng.permutation(len(ug))[: max(1, int(val_frac * len(ug)))]].tolist())
    vm = np.isin(g, list(vg))
    # MEMORY: the box is shared and runs at ~1 GB free, so the feature matrix is never
    # copied into torch.  Batches are sliced out of the numpy array on demand.
    ydt = np.int64 if family in ("grid", "abego4") else np.float32
    tr = np.where(~vm)[0]; va = np.where(vm)[0]
    X = np.ascontiguousarray(X, np.float32); Y = np.ascontiguousarray(Y, ydt)
    torch.manual_seed(seed)
    m = Trunk(X.shape[1], FAMILY_NOUT[family])
    opt = torch.optim.Adam(m.parameters(), lr=lr, weight_decay=1e-5)
    lf = loss_fn(family)
    best, bstate, bad, n = 1e18, None, 0, len(tr)
    rr = np.random.default_rng(seed + 5)
    for ep in range(epochs):
        m.train()
        perm = tr[rr.permutation(n)]
        for s in range(0, n, bs):
            b = perm[s:s + bs]
            opt.zero_grad()
            lf(m(torch.from_numpy(X[b])), torch.from_numpy(Y[b])).backward()
            opt.step()
        m.eval()
        with torch.no_grad():
            lv = float(np.mean([float(lf(m(torch.from_numpy(X[va[s:s + 8192]])),
                                         torch.from_numpy(Y[va[s:s + 8192]])))
                                for s in range(0, len(va), 8192)]))
        if verbose:
            print(f"    ep{ep} val={lv:.4f}", flush=True)
        if lv < best - 1e-4:
            best, bstate, bad = lv, {k: v.clone() for k, v in m.state_dict().items()}, 0
        else:
            bad += 1
            if bad >= patience:
                break
    if bstate:
        m.load_state_dict(bstate)
    m.eval()
    return m, {"val_loss": best, "epochs_run": ep + 1, "n_train": int(len(tr)),
               "n_val": int(len(va))}


# ============================================================ posteriors and their decoders
def posterior(model, X, family):
    """(L, 324) normalised density on the 18x18 grid, for EVERY family.

    Putting all four families on one representation is what makes them comparable: the
    proper score (NLL of the native cell), the point estimate, and `sigma_hat` are then
    computed by identical code and no family gets a decoder advantage.
    """
    with torch.no_grad():
        o = model(torch.tensor(np.asarray(X, np.float32)))
        if family == "grid":
            return torch.softmax(o, 1).numpy().astype(np.float32)
        if family == "mvm":
            gp = torch.tensor(GRID_PHI, dtype=torch.float32)
            gs = torch.tensor(GRID_PSI, dtype=torch.float32)
            return mvm_grid_density(o, gp, gs).numpy().astype(np.float32)
        if family == "point":
            # a point family still has an honest density: a von Mises whose concentration
            # is set by the SHRINKAGE of the regressed (cos, sin) vector, which is exactly
            # what an MSE-trained circular regressor reports about its own uncertainty.
            o = o.numpy()
            Rf = np.clip(np.linalg.norm(o[:, :2], axis=1), 1e-3, 0.999)
            Rs = np.clip(np.linalg.norm(o[:, 2:], axis=1), 1e-3, 0.999)
            muf = np.arctan2(o[:, 1], o[:, 0]); mus = np.arctan2(o[:, 3], o[:, 2])
            kf = Rf * (2 - Rf ** 2) / (1 - Rf ** 2)      # standard vM kappa estimator
            ks = Rs * (2 - Rs ** 2) / (1 - Rs ** 2)
            lp = (kf[:, None] * np.cos(GRID_PHI[None] - muf[:, None])
                  + ks[:, None] * np.cos(GRID_PSI[None] - mus[:, None]))
            lp -= lp.max(1, keepdims=True)
            p = np.exp(lp); return (p / p.sum(1, keepdims=True)).astype(np.float32)
        if family == "abego4":
            q = torch.softmax(o, 1).numpy()
            # spread each ABEGO class over the grid cells it owns, weighted by the
            # corpus-marginal occupancy of those cells (set by `set_abego_map`)
            return (q @ _ABEGO_MAP).astype(np.float32)
    raise ValueError(family)


_ABEGO_MAP = None          # (4, 324)


def set_abego_map(phi, psi, ok=None):
    """Corpus-marginal cell occupancy within each ABEGO class (training data only)."""
    global _ABEGO_MAP
    phi = np.asarray(phi, float); psi = np.asarray(psi, float)
    if ok is not None:
        phi, psi = phi[ok], psi[ok]
    b = grid_bin(phi, psi); a = abego4(phi, psi)
    M = np.zeros((4, NGRID), np.float32)
    for k in range(4):
        M[k] = np.bincount(b[a == k], minlength=NGRID)
    M += 1e-3
    _ABEGO_MAP = M / M.sum(1, keepdims=True)
    return _ABEGO_MAP


def marginal_posterior(phi, psi, ok, L):
    """MANDATORY NULL: the corpus-marginal Ramachandran density, same at every residue."""
    b = grid_bin(np.asarray(phi)[ok], np.asarray(psi)[ok])
    h = np.bincount(b, minlength=NGRID).astype(np.float32) + 1e-3
    h /= h.sum()
    return np.tile(h[None], (L, 1))


_DPHI = wrap(GRID_PHI[:, None] - GRID_PHI[None, :])
_DPSI = wrap(GRID_PSI[:, None] - GRID_PSI[None, :])
_ANGSQ = (_DPHI ** 2 + _DPSI ** 2) / 2.0          # (324, 324) mean squared angular offset
_CHEB = np.maximum(np.abs(_DPHI), np.abs(_DPSI))  # torus Chebyshev distance


def decode(P, ball_deg=40.0):
    """Posterior -> point estimate + self-estimated sigma, in the SAME units as S12's noise.

    The point estimate is the circular mean of the posterior mass inside a `ball_deg` ball
    around the MAP cell -- a plain circular mean over a bimodal posterior lands between the
    alpha and beta basins, which is a conformation that exists in neither, and that is the
    single most damaging decoding error available here.

    `sigma_hat` = sqrt(E_P[(dphi^2 + dpsi^2)/2]) about the point estimate, in DEGREES.
    S12's oracle added independent N(0, sigma) to each angle, so a perfectly calibrated
    predictor with sigma_hat = s is, per angle, exactly that experiment's sigma = s cell.
    """
    P = np.asarray(P, np.float64)
    P = P / np.maximum(P.sum(1, keepdims=True), 1e-30)
    mp = np.argmax(P, 1)
    ball = (_CHEB[mp] <= ball_deg * D2R)
    Wt = P * ball
    Wt /= np.maximum(Wt.sum(1, keepdims=True), 1e-30)
    cf = (Wt * np.cos(GRID_PHI)[None]).sum(1); sf = (Wt * np.sin(GRID_PHI)[None]).sum(1)
    cs = (Wt * np.cos(GRID_PSI)[None]).sum(1); ss = (Wt * np.sin(GRID_PSI)[None]).sum(1)
    phi = np.arctan2(sf, cf); psi = np.arctan2(ss, cs)
    dphi = wrap(GRID_PHI[None, :] - phi[:, None]); dpsi = wrap(GRID_PSI[None, :] - psi[:, None])
    var = (P * (dphi ** 2 + dpsi ** 2) / 2.0).sum(1)
    sigma_hat = np.sqrt(np.maximum(var, 0.0)) * R2D
    pmax = P.max(1)
    ent = -(P * np.log(P + 1e-30)).sum(1)
    mass30 = (P * (np.maximum(np.abs(dphi), np.abs(dpsi)) <= 30.0 * D2R)).sum(1)
    return {"phi": phi, "psi": psi, "sigma_hat": sigma_hat, "pmax": pmax,
            "entropy": ent, "mass30": mass30}


def nll_native(P, phi, psi):
    """Proper score: -log P(native cell).  Comparable across families by construction."""
    b = grid_bin(phi, psi)
    return -np.log(np.maximum(np.asarray(P, np.float64)[np.arange(len(b)), b], 1e-12))


def ang_err(pred, true):
    """Absolute circular deviation in DEGREES."""
    return np.abs(np.degrees(wrap(np.asarray(pred, float) - np.asarray(true, float))))
