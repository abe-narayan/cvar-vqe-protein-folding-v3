"""Pair-tensor distance predictor: the trRosetta/AlphaFold construction at peptide scale.

Why replace the per-pair MLP. `distogram.MLP` predicts every residue pair independently
from the concatenation of the two residues' ESM-2 embeddings and some scalars. That is the
crudest possible pair representation and it shows: correlation 0.744 with the true
distances, mean absolute error 2.04 A, in-band Spearman +0.425 against an oracle's +0.864.
It cannot know that d_ij, d_ik and d_kj are three sides of a triangle, that a beta ladder is
a diagonal band rather than a set of unrelated cells, or that a residue buried by one
contact must be near several others.

This module keeps a running (L, L, C) pair tensor and updates it with the three operations
that fix exactly those blindnesses:

* **triangle multiplicative update** -- `z_ij <- z_ij + sum_k a_ik * b_jk` (outgoing) and
  the incoming transpose. This is the operation that makes the representation of a pair
  depend on the two other sides of every triangle through it, and it is the single
  ingredient most responsible for AlphaFold's distograms being geometrically consistent.
  At L <= 26 it costs L^3 C, which is nothing.
* **axial attention** -- self-attention along rows and then columns, so a pair can read
  every other pair sharing a residue with it. This is what lets a ladder propagate.
* **transition** -- a per-pair MLP, the only thing the previous model had.

The output head is the same 17-bin distogram over the same bin edges, so everything
downstream (`distogram.Distogram`, the Bayes-risk score, the objective, the VQE) is
unchanged and the two model families are directly comparable on the same decoy bank.

Training is leave-fold-out over the same identity clusters, on the same peptides plus
fragments, so the comparison isolates the architecture.
"""
import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import distogram as dgm
import peptide_db as pdb
import priors

BASE = os.path.dirname(os.path.abspath(__file__))
import fragment_db as _fdb
MODEL_DIR = os.path.join(
    BASE, "pairnet_models_large" if _fdb.USE_LARGE else "pairnet_models")

MAXLEN = 26
NBINS = dgm.NBINS
BIN_EDGES = dgm.BIN_EDGES
CENTRES = dgm.CENTRES

#: Sequence-separation bins for the pair-feature encoding. Fine at short range, where the
#: separation nearly determines the distance, and coarse beyond.
SEP_BINS = np.array([1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 19, 24])
N_SEP = len(SEP_BINS) + 1


# ------------------------------------------------------------------ features
def residue_features(seq: str) -> np.ndarray:
    """``(L, D)`` per-residue inputs: ESM-2 embedding, physicochemistry, position."""
    n = len(seq)
    P = priors._prop(seq)
    pos = np.column_stack([np.arange(n) / max(1, n - 1),
                           np.minimum(np.arange(n), n - 1 - np.arange(n)) / max(1, n),
                           np.full(n, n / 26.0)])
    blocks = [P, pos]
    try:
        import esm_features
        if esm_features.available():
            blocks.insert(0, esm_features.embed(seq))
    except Exception:
        pass
    return np.hstack(blocks).astype(np.float32)


def pair_features(seq: str) -> np.ndarray:
    """``(L, L, D)`` per-pair inputs: separation encoding and the ESM-2 contact map."""
    n = len(seq)
    idx = np.arange(n)
    sep = np.abs(idx[:, None] - idx[None, :])
    onehot = np.zeros((n, n, N_SEP), np.float32)
    onehot[np.arange(n)[:, None].repeat(n, 1),
           np.arange(n)[None, :].repeat(n, 0),
           np.digitize(sep, SEP_BINS)] = 1.0
    extra = [onehot, (sep[..., None] / 26.0).astype(np.float32)]
    try:
        import esm_features
        if esm_features.available():
            con = esm_features.contacts(seq).astype(np.float32)
            extra.append(con[..., None])
            extra.append(np.log(np.clip(con, 1e-4, 1.0))[..., None] / 10.0)
    except Exception:
        pass
    return np.concatenate(extra, axis=-1)


def feature_dims(seq: str = "ACDEFGHIKL") -> Tuple[int, int]:
    return residue_features(seq).shape[1], pair_features(seq).shape[2]


# ------------------------------------------------------------------ model
def _build_module(d_res: int, d_pair: int, c: int, blocks: int, heads: int):
    import torch
    import torch.nn as nn

    class Triangle(nn.Module):
        """Triangle multiplicative update, outgoing and incoming."""

        def __init__(self, c, hidden=48):
            super().__init__()
            self.norm = nn.LayerNorm(c)
            self.a = nn.Linear(c, hidden)
            self.b = nn.Linear(c, hidden)
            self.ag = nn.Linear(c, hidden)
            self.bg = nn.Linear(c, hidden)
            self.og = nn.Linear(c, c)
            self.out = nn.Linear(hidden, c)
            self.onorm = nn.LayerNorm(hidden)

        def forward(self, z, mask):
            # z (B, L, L, C), mask (B, L, L)
            zn = self.norm(z)
            m = mask.unsqueeze(-1)
            a = torch.sigmoid(self.ag(zn)) * self.a(zn) * m
            b = torch.sigmoid(self.bg(zn)) * self.b(zn) * m
            out = torch.einsum("bikc,bjkc->bijc", a, b)          # outgoing
            out = out + torch.einsum("bkic,bkjc->bijc", a, b)    # incoming
            return torch.sigmoid(self.og(zn)) * self.out(self.onorm(out))

    class Axial(nn.Module):
        """Row-then-column self-attention over the pair tensor."""

        def __init__(self, c, heads):
            super().__init__()
            self.h = heads
            self.dh = c // heads
            self.norm = nn.LayerNorm(c)
            self.qkv = nn.Linear(c, 3 * c, bias=False)
            self.proj = nn.Linear(c, c)

        def _attend(self, x, bias):
            B, R, L, C = x.shape
            q, k, v = self.qkv(x).chunk(3, dim=-1)
            shape = (B, R, L, self.h, self.dh)
            q, k, v = (t.reshape(shape).permute(0, 1, 3, 2, 4) for t in (q, k, v))
            att = (q @ k.transpose(-1, -2)) / math.sqrt(self.dh) + bias
            att = att.softmax(-1)
            o = (att @ v).permute(0, 1, 3, 2, 4).reshape(B, R, L, C)
            return self.proj(o)

        def forward(self, z, mask):
            zn = self.norm(z)
            bias = (1.0 - mask[:, :1, :]).unsqueeze(1).unsqueeze(1) * -1e4
            row = self._attend(zn, bias)
            col = self._attend(zn.transpose(1, 2), bias).transpose(1, 2)
            return row + col

    class PairNet(nn.Module):
        def __init__(self):
            super().__init__()
            self.res_a = nn.Linear(d_res, c)
            self.res_b = nn.Linear(d_res, c)
            self.pair = nn.Linear(d_pair, c)
            self.blocks = nn.ModuleList()
            for _ in range(blocks):
                self.blocks.append(nn.ModuleList([
                    Triangle(c), Axial(c, heads),
                    nn.Sequential(nn.LayerNorm(c), nn.Linear(c, 2 * c), nn.GELU(),
                                  nn.Linear(2 * c, c))]))
            self.head = nn.Sequential(nn.LayerNorm(c), nn.Linear(c, NBINS))

        def forward(self, r, p, mask):
            z = self.res_a(r)[:, :, None, :] + self.res_b(r)[:, None, :, :] + self.pair(p)
            m = mask.unsqueeze(-1)
            z = z * m
            for tri, axi, tra in self.blocks:
                z = (z + tri(z, mask)) * m
                z = (z + axi(z, mask)) * m
                z = (z + tra(z)) * m
            logit = self.head(z)
            return 0.5 * (logit + logit.transpose(1, 2))          # symmetric by design

    return PairNet()


class PairNetModel:
    """Trainable wrapper: padding, masking, normalisation, soft-binned loss."""

    def __init__(self, d_res: int, d_pair: int, c: int = 64, blocks: int = 4,
                 heads: int = 4, seed: int = 0):
        import torch
        torch.manual_seed(seed)
        self.torch = torch
        self.d_res, self.d_pair = d_res, d_pair
        self.cfg = dict(c=c, blocks=blocks, heads=heads)
        self.net = _build_module(d_res, d_pair, c, blocks, heads)
        self.rmu = np.zeros(d_res, np.float32)
        self.rsd = np.ones(d_res, np.float32)

    # -- data -----------------------------------------------------------
    @staticmethod
    def encode(entries: Sequence[pdb.Peptide], with_target: bool = True):
        R, P, M, Y = [], [], [], []
        for e in entries:
            n = e.n
            r = np.zeros((MAXLEN, residue_features(e.seq).shape[1]), np.float32)
            rf = residue_features(e.seq)
            r[:n] = rf
            pf = pair_features(e.seq)
            p = np.zeros((MAXLEN, MAXLEN, pf.shape[2]), np.float16)
            p[:n, :n] = pf
            m = np.zeros((MAXLEN, MAXLEN), np.float32)
            m[:n, :n] = 1.0
            R.append(r)
            P.append(p)
            M.append(m)
            if with_target:
                y = np.zeros((MAXLEN, MAXLEN), np.int16)
                d = np.linalg.norm(e.ca[:, None] - e.ca[None], axis=-1)
                y[:n, :n] = np.digitize(d, BIN_EDGES)
                Y.append(y)
        out = [np.stack(R), np.stack(P), np.stack(M)]
        if with_target:
            out.append(np.stack(Y))
        return out

    # -- training -------------------------------------------------------
    def fit(self, entries: Sequence[pdb.Peptide], epochs: int = 60, batch: int = 24,
            lr: float = 1.2e-3, smooth: float = 0.6, verbose: bool = True):
        t = self.torch
        R, P, M, Y = self.encode(entries)
        # Normalise over REAL residues only; padding rows are zero and would otherwise
        # drag the mean toward zero in proportion to how much padding a batch carries.
        real = M[:, :, 0] > 0
        used = R[real]
        self.rmu = used.mean(0)
        self.rsd = used.std(0) + 1e-3
        R = (R - self.rmu) / self.rsd
        idx = np.arange(NBINS)
        soft = np.exp(-((idx[None, :] - idx[:, None]) ** 2) / (2 * smooth ** 2))
        soft = (soft / soft.sum(1, keepdims=True)).astype(np.float32)
        soft_t = t.as_tensor(soft)
        # Kept as numpy; each batch is cast on the way in, so the float32 copy is one
        # batch rather than the whole set.
        Rt, Mt = t.as_tensor(R), t.as_tensor(M)
        Pt, Yt = P, Y
        opt = t.optim.AdamW(self.net.parameters(), lr=lr, weight_decay=1e-4)
        per_epoch = int(np.ceil(len(R) / batch))
        sched = t.optim.lr_scheduler.OneCycleLR(
            opt, max_lr=lr, total_steps=epochs * per_epoch + 1)
        # Length bucketing: batches are formed from chains of similar length and the
        # padded tensors cropped to the batch maximum. The triangle update costs L^3, so
        # padding a batch of 9-mers out to 26 wastes 24x on that term; sorting by length
        # and cropping is a ~2x end-to-end speedup with no change to what is computed
        # (padded rows are masked out either way).
        lengths = np.asarray([e.n for e in entries])
        order = np.argsort(lengths, kind="stable")
        groups = [order[i:i + batch] for i in range(0, len(order), batch)]
        crops = [int(lengths[g].max()) for g in groups]
        g = t.Generator().manual_seed(0)
        self.net.train()
        for ep in range(epochs):
            gperm = t.randperm(len(groups), generator=g).tolist()
            tot = w = 0.0
            for gi in gperm:
                b = t.as_tensor(groups[gi])
                L = crops[gi]
                bi = b.numpy()
                logit = self.net(Rt[b][:, :L],
                                 t.as_tensor(Pt[bi][:, :L, :L].astype(np.float32)),
                                 Mt[b][:, :L, :L])
                target = soft_t[t.as_tensor(Yt[bi][:, :L, :L].astype(np.int64))]
                lp = t.log_softmax(logit, -1)
                loss = -(target * lp).sum(-1)
                m = Mt[b][:, :L, :L]
                loss = (loss * m).sum() / m.sum().clamp(min=1.0)
                opt.zero_grad()
                loss.backward()
                t.nn.utils.clip_grad_norm_(self.net.parameters(), 1.0)
                opt.step()
                sched.step()
                tot += float(loss.detach()) * float(m.sum())
                w += float(m.sum())
            if verbose and ep % 10 == 0:
                print(f"    pairnet epoch {ep} loss {tot / max(w, 1):.4f}", flush=True)
        return self

    # -- inference ------------------------------------------------------
    def predict_map(self, seq: str) -> np.ndarray:
        """``(n, n, NBINS)`` distance-bin probabilities for one sequence."""
        t = self.torch
        n = len(seq)
        r = np.zeros((1, MAXLEN, self.d_res), np.float32)
        r[0, :n] = residue_features(seq)
        r = (r - self.rmu) / self.rsd
        pf = pair_features(seq)
        p = np.zeros((1, MAXLEN, MAXLEN, self.d_pair), np.float32)
        p[0, :n, :n] = pf
        m = np.zeros((1, MAXLEN, MAXLEN), np.float32)
        m[0, :n, :n] = 1.0
        self.net.eval()
        with t.no_grad():
            logit = self.net(t.as_tensor(r), t.as_tensor(p), t.as_tensor(m))
            prob = t.softmax(logit, -1).numpy()[0, :n, :n]
        return prob

    def save(self, path: str):
        self.torch.save({"sd": self.net.state_dict(), "rmu": self.rmu, "rsd": self.rsd,
                         "cfg": self.cfg, "d": (self.d_res, self.d_pair)}, path)

    def load(self, path: str):
        z = self.torch.load(path, weights_only=False)
        self.net.load_state_dict(z["sd"])
        self.rmu, self.rsd = z["rmu"], z["rsd"]
        return self


# ------------------------------------------------------------------ folds
def train_fold(fold: int, n_folds: int = 5, fragments: bool = True,
               epochs: int = 60, c: int = 64, blocks: int = 4, seed: int = 0,
               verbose: bool = True) -> PairNetModel:
    os.makedirs(MODEL_DIR, exist_ok=True)
    path = os.path.join(MODEL_DIR, f"fold{fold}_c{c}b{blocks}_s{seed}.pt")
    d_res, d_pair = feature_dims()
    m = PairNetModel(d_res, d_pair, c=c, blocks=blocks, seed=seed)
    if os.path.exists(path):
        return m.load(path)
    entries = [p for p in pdb.load() if pdb.folds(n_folds)[p.seq] != fold]
    if fragments:
        entries = entries + dgm._fold_fragments(fold, n_folds)
    if verbose:
        print(f"  pairnet fold {fold}: {len(entries)} chains, "
              f"d_res={d_res} d_pair={d_pair} c={c} blocks={blocks}", flush=True)
    m.fit(entries, epochs=epochs, verbose=verbose)
    m.save(path)
    return m


def distogram_for(sequence: str, fold: Optional[int] = None, n_folds: int = 5,
                  model: Optional[PairNetModel] = None, **kw) -> dgm.Distogram:
    """A `distogram.Distogram` whose probabilities come from the pair network."""
    if model is None:
        if fold is None:
            fold = pdb.folds(n_folds).get(sequence, 0)
        model = train_fold(fold, n_folds, verbose=False, **kw)
    n = len(sequence)
    full = model.predict_map(sequence)
    i, j = np.triu_indices(n, k=2)
    prob = np.maximum(full[i, j], 1e-6)
    prob = prob / prob.sum(1, keepdims=True)
    return dgm.Distogram(sequence, prob.astype(np.float64), i, j)
