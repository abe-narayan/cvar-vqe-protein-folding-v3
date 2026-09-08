"""EXPERIMENT 2 -- the learned SET DECODER.

A permutation-invariant model that reads the FULL (m x npairs x F) candidate-vs-objective
deviation map plus the (m x m) candidate-vs-candidate RMSD matrix and emits weights over
the candidate set.  Trained leave-fold-out with an objective that is differentiable through
the coordinate average (NOT through the projection):

    L = CA-RMSD( sum_c w_c Asup_c , native )        [rotation detached: envelope theorem]

Comparison heads: `listwise` (CE against softmax(-rr/T)) and `regress` (MSE on rr, the class
of objective every previous scalar-feature ranker used).

Discipline: leave-fold-out (5 pinned folds), inner-fold early stopping, learning curve over
training-set size, a feature-permutation null, and an oracle-feature positive control.

Usage:  python -m s12.agg_decoder <config-name> [...]
"""
from __future__ import annotations
import os, sys, json, time, math
import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
import torch
import torch.nn as nn
import torch.nn.functional as Fn

torch.set_num_threads(int(os.environ.get("TORCH_THREADS", "1")))   # 2 threads thrash on this box

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from s12 import instrument as I
from s12 import agg_features as AF

OUT = os.path.join(ROOT, "s12", "cache", "agg_dec")
os.makedirs(OUT, exist_ok=True)


# ------------------------------------------------------------------ data
def load_all(source="top75", abs_feat=False):
    tg = I.targets()
    data = []
    for t in tg:
        o = AF.load(t["pdb"], source)
        data.append({"pdb": t["pdb"], "fold": t["fold"], "n": t["n"],
                     "X": torch.tensor(o["X"], dtype=torch.float32),
                     "G": torch.tensor(o["G"], dtype=torch.float32),
                     "P": torch.tensor(o["P"], dtype=torch.float32),
                     "A": torch.tensor(o["Asup"], dtype=torch.float32),
                     "nat": torch.tensor(o["nat"], dtype=torch.float32),
                     "rr": torch.tensor(o["rr"], dtype=torch.float32),
                     "sc": torch.tensor(o["sc"], dtype=torch.float32)})
    if abs_feat:
        data = [abs_signed(d) for d in data]
    return data


def norm_stats(data, idx):
    xs = torch.cat([data[i]["X"].reshape(-1, data[i]["X"].shape[-1]) for i in idx], 0)
    gs = torch.cat([data[i]["G"] for i in idx], 0)
    return (xs.mean(0), xs.std(0) + 1e-6, gs.mean(0), gs.std(0) + 1e-6)


# ------------------------------------------------------------------ differentiable RMSD
def krmsd(C, T):
    """CA-RMSD(C, T) with the optimal rotation detached (envelope theorem => exact grad)."""
    Cc = C - C.mean(0, keepdim=True)
    Tc = T - T.mean(0, keepdim=True)
    with torch.no_grad():
        H = Cc.t() @ Tc
        U, S, Vh = torch.linalg.svd(H)
        d = torch.sign(torch.det(Vh.t() @ U.t()))
        D = torch.diag(torch.tensor([1.0, 1.0, 1.0]))
        D[2, 2] = d
        R = Vh.t() @ D @ U.t()
    return torch.sqrt((((Cc @ R.t()) - Tc) ** 2).sum() / C.shape[0] + 1e-9)


# ------------------------------------------------------------------ model
class SAB(nn.Module):
    def __init__(self, H, heads):
        super().__init__()
        self.h = heads
        self.q = nn.Linear(H, H); self.k = nn.Linear(H, H); self.v = nn.Linear(H, H)
        self.o = nn.Linear(H, H)
        self.ln1 = nn.LayerNorm(H); self.ln2 = nn.LayerNorm(H)
        self.ff = nn.Sequential(nn.Linear(H, 2 * H), nn.GELU(), nn.Linear(2 * H, H))
        self.alpha = nn.Parameter(torch.zeros(heads))

    def forward(self, x, P):
        m, H = x.shape
        dh = H // self.h
        q = self.q(x).view(m, self.h, dh).transpose(0, 1)
        k = self.k(x).view(m, self.h, dh).transpose(0, 1)
        v = self.v(x).view(m, self.h, dh).transpose(0, 1)
        att = (q @ k.transpose(1, 2)) / math.sqrt(dh)
        att = att - Fn.softplus(self.alpha).view(-1, 1, 1) * P.unsqueeze(0)
        a = torch.softmax(att, -1)
        y = (a @ v).transpose(0, 1).reshape(m, H)
        x = self.ln1(x + self.o(y))
        return self.ln2(x + self.ff(x))


class SetDecoder(nn.Module):
    def __init__(self, F=AF.F, GF=AF.GF, H=32, heads=4, layers=2, drop=0.1):
        super().__init__()
        self.pair = nn.Sequential(nn.Linear(F, H), nn.GELU(), nn.Linear(H, H), nn.GELU())
        self.proj = nn.Sequential(nn.Linear(3 * H + GF, H), nn.GELU(), nn.Dropout(drop))
        self.sab = nn.ModuleList([SAB(H, heads) for _ in range(layers)])
        self.head = nn.Sequential(nn.Linear(H, H), nn.GELU(), nn.Linear(H, 1))
        self.logtau = nn.Parameter(torch.zeros(1))

    def forward(self, X, G, P):
        h = self.pair(X)                                   # (m, npairs, H)
        z = torch.cat([h.mean(1), h.max(1).values, h.std(1), G], -1)
        z = self.proj(z)
        for s in self.sab:
            z = s(z, P)
        return self.head(z).squeeze(-1)                    # (m,) logits


# ------------------------------------------------------------------ train / eval
def prep(d, st, permute=None, oracle_feat=False, zero_feat=False, gen=None):
    xm, xs, gm, gs = st
    X = (d["X"] - xm) / xs
    G = (d["G"] - gm) / gs
    if oracle_feat:
        rr = d["rr"]
        G = torch.cat([G, ((rr - rr.mean()) / (rr.std() + 1e-6)).unsqueeze(-1)], -1)
    if permute is not None:
        X = X[permute]; G = G[permute]
    if zero_feat:
        X = torch.zeros_like(X); G = torch.zeros_like(G)
    return X, G, d["P"]


def abs_signed(d):
    """Magnitude-only ablation: |.| on exactly the SIGNED per-pair deviation columns.
    Isolates whether the SIGN of the deviation map carries the information."""
    X = d["X"].clone()
    X[:, :, AF.SIGNED_PAIR_COLS] = X[:, :, AF.SIGNED_PAIR_COLS].abs()
    d = dict(d); d["X"] = X
    return d


def emit(d, w):
    return (d["A"] * w.view(-1, 1, 1)).sum(0)


def run_fold(data, test_fold, cfg, seed=0, ntrain=None, log=None):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    folds = np.array([d["fold"] for d in data])
    te = np.where(folds == test_fold)[0]
    tr_all = np.where(folds != test_fold)[0]
    if cfg.get("no_fail18_train"):
        fset = set(I.FAIL18)
        tr_all = np.array([i for i in tr_all if data[i]["pdb"] not in fset])
    inner = sorted(set(folds[tr_all].tolist()))
    preds = {}
    models = 0
    acc = {int(i): [] for i in te}
    for gval in inner[: cfg.get("n_inner", 1)]:
        if cfg.get("use_all"):
            tr = np.array(tr_all); va = np.array(tr_all)
        else:
            tr = np.array([i for i in tr_all if folds[i] != gval])
            va = np.array([i for i in tr_all if folds[i] == gval])
        if ntrain is not None and ntrain < len(tr):
            tr = rng.choice(tr, ntrain, replace=False)
        st = norm_stats(data, tr)
        gf = AF.GF + (1 if cfg.get("oracle_feat") else 0)
        model = SetDecoder(F=AF.F, GF=gf, H=cfg.get("H", 32), heads=cfg.get("heads", 4),
                           layers=cfg.get("layers", 2), drop=cfg.get("drop", 0.1))
        opt = torch.optim.AdamW(model.parameters(), lr=cfg.get("lr", 3e-4), weight_decay=cfg.get("wd", 1e-4))
        perm = {}
        if cfg.get("null"):
            for i in range(len(data)):
                perm[i] = torch.tensor(rng.permutation(len(data[i]["rr"])))
        if cfg.get("null_A"):
            for i in range(len(data)):
                q = torch.tensor(rng.permutation(len(data[i]["rr"])))
                data[i]["A"] = data[i]["A"][q]
        best, best_state, bad = 1e9, None, 0
        for ep in range(cfg.get("epochs", 80)):
            model.train()
            order = rng.permutation(len(tr))
            opt.zero_grad()
            for bi, k in enumerate(order):
                d = data[tr[k]]
                X, G, P = prep(d, st, perm.get(int(tr[k])), cfg.get("oracle_feat"), cfg.get("zero_feat"))
                s = model(X, G, P)
                loss = _loss(s, d, cfg)
                (loss / cfg.get("accum", 8)).backward()
                if (bi + 1) % cfg.get("accum", 8) == 0 or bi == len(order) - 1:
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                    opt.step(); opt.zero_grad()
            model.eval()
            with torch.no_grad():
                v = float(np.mean([float(krmsd(emit(data[i], _w(model(*prep(data[i], st, perm.get(int(i)), cfg.get("oracle_feat"), cfg.get("zero_feat"))), cfg)),
                                              data[i]["nat"])) for i in va]))
            if v < best - 1e-4:
                best, bad = v, 0
                best_state = {k2: v2.detach().clone() for k2, v2 in model.state_dict().items()}
            else:
                bad += 1
                if bad >= cfg.get("patience", 15):
                    break
            if log is not None and ep % 10 == 0:
                print(f"  f{test_fold} g{gval} ep{ep} val {v:.4f} best {best:.4f}", flush=True)
        model.load_state_dict(best_state)
        model.eval()
        with torch.no_grad():
            for i in te:
                s = model(*prep(data[i], st, perm.get(int(i)), cfg.get("oracle_feat"), cfg.get("zero_feat")))
                acc[int(i)].append(_w(s, cfg).numpy())
        models += 1
    for i in te:
        preds[int(i)] = np.mean(acc[int(i)], 0)
    return preds


def _w(s, cfg):
    if cfg.get("head", "weight") == "regress":
        # predicted rr -> softmax(-z/T)
        z = (s - s.mean()) / (s.std() + 1e-6)
        return torch.softmax(-z / cfg.get("wtemp", 0.5), 0)
    return torch.softmax(s / torch.tensor(1.0), 0)


def _loss(s, d, cfg):
    head = cfg.get("head", "weight")
    if head == "weight":
        w = torch.softmax(s, 0)
        return krmsd(emit(d, w), d["nat"])
    if head == "listwise":
        rr = d["rr"]
        tgt = torch.softmax(-(rr - rr.min()) / cfg.get("ltemp", 0.3), 0)
        return -(tgt * torch.log_softmax(s, 0)).sum()
    rr = d["rr"]
    return Fn.mse_loss(s, (rr - rr.mean()) / (rr.std() + 1e-6))


def evaluate(data, preds):
    out = {}
    for i, w in preds.items():
        d = data[i]
        C = (d["A"].numpy() * np.asarray(w)[:, None, None]).sum(0)
        out[d["pdb"]] = {"raw": float(I.ca_rmsd(C, d["nat"].numpy())),
                         "w": np.asarray(w, float).tolist(),
                         "ess": float(1.0 / (np.asarray(w, float) ** 2).sum()),
                         "rho": float(np.corrcoef(np.asarray(w, float), d["rr"].numpy())[0, 1])}
    return out


def main(name, cfg, ntrain=None, seed=0, verbose=True):
    data = load_all(cfg.get("source", "top75"), abs_feat=cfg.get("abs_feat", False))
    preds = {}
    t0 = time.time()
    for f in sorted(set(d["fold"] for d in data)):
        preds.update(run_fold(data, f, cfg, seed=seed, ntrain=ntrain, log=verbose or None))
        if verbose:
            print(f"fold {f} done {time.time()-t0:.0f}s", flush=True)
    ev = evaluate(data, preds)
    np.savez_compressed(os.path.join(OUT, f"{name}.npz"),
                        **{p: np.asarray(v["w"]) for p, v in ev.items()})
    res = {"cfg": cfg, "ntrain": ntrain, "seed": seed,
           "raw": {p: v["raw"] for p, v in ev.items()},
           "ess": {p: v["ess"] for p, v in ev.items()},
           "rho": {p: v["rho"] for p, v in ev.items()},
           "mean_raw": float(np.mean([v["raw"] for v in ev.values()])),
           "mean_rho": float(np.nanmean([v["rho"] for v in ev.values()])),
           "secs": time.time() - t0}
    I.write(f"agg_dec_{name}", res)
    print(name, "mean_raw", round(res["mean_raw"], 4), "mean_rho", round(res["mean_rho"], 3),
          "secs", int(res["secs"]), flush=True)
    return res


CFGS = {
    "v2":        {"head": "weight", "no_fail18_train": True},
    "v2_abs":    {"head": "weight", "no_fail18_train": True, "abs_feat": True},
    "v2_null":   {"head": "weight", "no_fail18_train": True, "null": True},
    "v2_oracle": {"head": "weight", "no_fail18_train": True, "oracle_feat": True},
    "v2_101":    {"head": "weight", "no_fail18_train": True, "use_all": True, "epochs": 40, "patience": 99},
    "v2_hi":     {"head": "weight", "no_fail18_train": True, "H": 64, "layers": 3, "epochs": 150, "patience": 30, "lr": 1e-3},
    "weight101": {"head": "weight", "use_all": True, "epochs": 40, "patience": 99},
    "weight_s1": {"head": "weight"},
    "oracle101": {"head": "weight", "oracle_feat": True, "use_all": True, "epochs": 40, "patience": 99},
    "weight_hi": {"head": "weight", "H": 64, "layers": 3, "epochs": 150, "patience": 30, "lr": 1e-3},
    "weight": {"head": "weight"},
    "listwise": {"head": "listwise"},
    "regress": {"head": "regress"},
    "null": {"head": "weight", "null": True},
    "oracle": {"head": "weight", "oracle_feat": True},
    "weight_deep": {"head": "weight", "H": 64, "layers": 3},
    "weight_noattn": {"head": "weight", "layers": 0},
    "null_A": {"head": "weight", "null_A": True},
    "weight_p_only": {"head": "weight", "zero_feat": True},
}

if __name__ == "__main__":
    args = sys.argv[1:]
    nm = args[0]
    cfg = dict(CFGS[nm])
    nt = int(args[1]) if len(args) > 1 and args[1] != "-" else None
    sd = int(args[2]) if len(args) > 2 else 0
    tag = nm + (f"_n{nt}" if nt else "") + (f"_s{sd}" if sd else "")
    main(tag, cfg, ntrain=nt, seed=sd)
