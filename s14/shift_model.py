"""SHIFT agent -- STEPS 3-5.  A REAL shift -> torsion predictor, and the asymmetry test
that decides whether the experimental channel carries anything the sequence does not.

This is not a simulation.  It trains a TALOS-class model on (backbone chemical shift ->
phi/psi) pairs harvested from BMRB depositions with matched PDB coordinates (s14/
shift_corpus.py), all of them outside the 126-target instrument and containment-screened
against every target sequence, and applies it to the MEASURED shifts of the targets where
STEP 1 found them.

The representation is deliberately the same 18x18 = 324-cell Ramachandran grid that
`s13/tors_common.py` already uses -- which is also, verbatim, TALOS-N's own output
representation (Shen & Bax 2013 eq. 1-2: "the 360 x 360 Ramachandran map is binned into
18 x 18 square boxes, or voxels ... a backbone phi/psi distribution code with 324 states").
So the shift arm, the sequence-only arm and TALOS-N all speak one language and the decoder
is shared code.

ARMS, each labelled with exactly what information it had:

  E-SHIFT    measured BMRB shifts of the target        -- NMR-RESTRAINED, own column
  P-SHIFT    shifts predicted from sequence alone      -- sequence-only in disguise
  REF-SHIFT  all secondary shifts set to zero          -- NULL: zero experimental content
  MASK-ONLY  real assignment pattern, shift values zeroed -- NULL: does the PATTERN leak?
  SEQ-ONLY   s13's best sequence-only arm (a_pepPos)   -- the incumbent sequence channel

    python -m s14.shift_model
"""
from __future__ import annotations

import os
import sys
import json
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import torch                                            # noqa: E402
import torch.nn as nn                                   # noqa: E402

from s12 import instrument as I                         # noqa: E402
from s13 import tors_common as T                        # noqa: E402
from s13 import tors_eval as EV                         # noqa: E402
from s13.tors_train import post_path                    # noqa: E402
from s14.shift_bmrb import BACKBONE, RESULTS            # noqa: E402
from s14 import shift_bmrb as B                         # noqa: E402
from s14.shift_avail import position_weights            # noqa: E402

torch.set_num_threads(2)

AA = "ACDEFGHIKLMNPQRSTVWY"
AAI = {a: i for i, a in enumerate(AA)}
WIN = (-1, 0, 1)
NF_SH = len(WIN) * len(BACKBONE) * 2                    # value + present flag
NF_AA = len(WIN) * (len(AA) + 1)
NFEAT = NF_SH + NF_AA
D2R = np.pi / 180.0
R2D = 180.0 / np.pi
SEED = 0


# ------------------------------------------------------------------ corpus -> features
def load_corpus():
    with open(os.path.join(RESULTS, "shift_corpus.json")) as fh:
        return json.load(fh)


def reference(rows):
    """Per (residue type, nucleus) reference shift, as the corpus MEDIAN, plus a per-nucleus
    scale.  This is a corpus reference, NOT a published random-coil table; it is a pure
    function of residue type, so subtracting it removes no experimental content."""
    acc = {}
    for r in rows:
        for a, v in r["sh"].items():
            acc.setdefault((r["aa"], a), []).append(v)
    ref = {k: float(np.median(v)) for k, v in acc.items() if len(v) >= 8}
    glob = {}
    for a in BACKBONE:
        v = [x for (aa, nu), lst in acc.items() if nu == a for x in lst]
        glob[a] = float(np.median(v)) if v else 0.0
    sc = {}
    for a in BACKBONE:
        d = [v - ref.get((aa, nu), glob[a])
             for (aa, nu), lst in acc.items() if nu == a for v in lst]
        sc[a] = float(np.std(d)) if d else 1.0
        if sc[a] < 1e-6:
            sc[a] = 1.0
    return ref, glob, sc


def feat_row(aas, shs, i, ref, glob, sc, zero_values=False, force_present=False):
    """aas: list of one-letter codes (or '-'), shs: list of {nucleus: value} per position."""
    f = np.zeros(NFEAT, np.float32)
    k = 0
    for w in WIN:
        j = i + w
        ok = 0 <= j < len(aas)
        aa = aas[j] if ok else "-"
        sh = shs[j] if ok else {}
        for a in BACKBONE:
            v = sh.get(a)
            if v is None:
                f[k] = 0.0
                f[k + 1] = 1.0 if force_present else 0.0
            else:
                r = ref.get((aa, a), glob[a])
                f[k] = 0.0 if zero_values else float((v - r) / sc[a])
                f[k + 1] = 1.0
            k += 2
    for w in WIN:
        j = i + w
        aa = aas[j] if 0 <= j < len(aas) else "-"
        idx = AAI.get(aa, len(AA))
        f[k + idx] = 1.0
        k += len(AA) + 1
    return f


def corpus_xy(rows, ref, glob, sc):
    X = np.zeros((len(rows), NFEAT), np.float32)
    phi = np.zeros(len(rows)); psi = np.zeros(len(rows))
    for n, r in enumerate(rows):
        aas = [r["prev"], r["aa"], r["next"]]
        shs = [r["sh_prev"], r["sh"], r["sh_next"]]
        X[n] = feat_row(aas, shs, 1, ref, glob, sc)
        phi[n] = r["phi"]; psi[n] = r["psi"]
    return X, phi, psi


# ------------------------------------------------------------------------------ models
class MLP(nn.Module):
    def __init__(self, nin, nout, hid=384, p=0.2):
        super().__init__()
        self.f = nn.Sequential(
            nn.Linear(nin, hid), nn.ReLU(), nn.Dropout(p),
            nn.Linear(hid, hid), nn.ReLU(), nn.Dropout(p),
            nn.Linear(hid, nout))

    def forward(self, x):
        return self.f(x)


def fit(X, y, nout, task="cls", epochs=40, bs=512, lr=1e-3, seed=SEED, val=0.1,
        groups=None, verbose=True, tag=""):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    n = len(X)
    if groups is None:
        idx = rng.permutation(n)
    else:                                   # group-held-out validation (by protein)
        gs = np.array(groups)
        ug = rng.permutation(np.unique(gs))
        nv = max(1, int(val * len(ug)))
        vg = set(ug[:nv].tolist())
        vi = np.flatnonzero(np.isin(gs, list(vg)))
        ti = np.flatnonzero(~np.isin(gs, list(vg)))
        idx = np.concatenate([vi, ti])
        nval = len(vi)
    if groups is None:
        nval = max(1, int(val * n))
    Xv = torch.tensor(X[idx[:nval]]); Xt = torch.tensor(X[idx[nval:]])
    if task == "cls":
        yv = torch.tensor(y[idx[:nval]]); yt = torch.tensor(y[idx[nval:]])
        lossf = nn.CrossEntropyLoss()
    else:
        yv = torch.tensor(y[idx[:nval]]); yt = torch.tensor(y[idx[nval:]])
        lossf = None
    m = MLP(X.shape[1], nout)
    opt = torch.optim.Adam(m.parameters(), lr=lr)
    best, bstate, bad = 1e18, None, 0
    for ep in range(epochs):
        m.train()
        perm = torch.randperm(len(Xt))
        for s in range(0, len(Xt), bs):
            b = perm[s:s + bs]
            opt.zero_grad()
            o = m(Xt[b])
            if task == "cls":
                L = lossf(o, yt[b])
            else:
                msk = ~torch.isnan(yt[b])
                L = (((o - torch.nan_to_num(yt[b])) ** 2) * msk).sum() / msk.sum().clamp(min=1)
            L.backward(); opt.step()
        m.eval()
        with torch.no_grad():
            o = m(Xv)
            if task == "cls":
                v = float(lossf(o, yv))
            else:
                msk = ~torch.isnan(yv)
                v = float((((o - torch.nan_to_num(yv)) ** 2) * msk).sum()
                          / msk.sum().clamp(min=1))
        if v < best - 1e-4:
            best, bstate, bad = v, {k: t.clone() for k, t in m.state_dict().items()}, 0
        else:
            bad += 1
            if bad >= 6:
                break
    if bstate:
        m.load_state_dict(bstate)
    if verbose:
        print("    fit {:10s} n={:6d} nout={:4d} best_val={:.4f} epochs={}".format(
            tag, n, nout, best, ep + 1), flush=True)
    m.eval()
    return m


def posterior(m, X):
    with torch.no_grad():
        o = m(torch.tensor(X))
        return torch.softmax(o, 1).numpy().astype(np.float64)


# -------------------------------------------------------------------- target-side input
def target_shifts(pdbid, A):
    """-> (aas, shs, gate_mask).  shs[i] = {nucleus: value} of MEASURED shifts."""
    a = A[pdbid]
    b = a["best"]
    seq = a["seq"]
    aas = list(seq)
    if b is None:
        return aas, [{} for _ in aas], np.zeros(len(aas), bool)
    shs = []
    for i in range(len(aas)):
        d = {}
        for nu in BACKBONE:
            if b["per_nuc"][nu][i]:
                d[nu] = None                     # placeholder, filled below
        shs.append(d)
    return aas, shs, np.array(b["talos_mask"], bool)


def target_shift_values(pdbid, A):
    """Re-read the actual ppm values for the target's window from the cached BMRB entry."""
    a = A[pdbid]
    b = a["best"]
    seq = a["seq"]
    if b is None:
        return [{} for _ in seq]
    rows = B.entry_shifts(b["bmrb"])
    ents = B.entry_entities(b["bmrb"])
    ent = [e for e in ents if e["id"] == b["entity"]][0]
    loc = B.locate(seq, ent["seq"])
    have = {}
    for r in rows:
        if str(r["entity"]) != str(b["entity"]):
            continue
        nu = B.ATOM_ALIAS.get(r["atom"])
        if nu is None:
            continue
        have.setdefault(r["seq_id"], {})[nu] = r["val"]
    return [dict(have.get(s, {})) if s is not None else {} for s in loc["map"]]


# ------------------------------------------------------------------------------ scoring
def sigma_eff(pred_phi, pred_psi, nphi, npsi, det_mask):
    dp = T.wrap(pred_phi - nphi)[det_mask[0]]
    ds = T.wrap(pred_psi - npsi)[det_mask[1]]
    e = np.concatenate([dp, ds])
    return float(np.sqrt((e ** 2).mean()) * R2D)


def lag1(errs):
    """Along-chain lag-1 autocorrelation of the signed torsion error -- the coherence axis
    the coordinator's s14/coherence.py parameterises the restraint surface by."""
    a, b = [], []
    for e in errs:
        if len(e) < 3:
            continue
        a.append(e[:-1]); b.append(e[1:])
    if not a:
        return float("nan")
    a = np.concatenate(a); b = np.concatenate(b)
    if a.std() < 1e-9 or b.std() < 1e-9:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def main():
    t0 = time.time()
    T.free_ok(1.0, tag="shift_model")
    rows = load_corpus()
    prots = sorted({r["pdb"] for r in rows})
    print("corpus: {} proteins, {} residues".format(len(prots), len(rows)))

    with open(os.path.join(RESULTS, "shift_availability.json")) as fh:
        A = json.load(fh)
    with open(os.path.join(RESULTS, "shift_avail_summary.json")) as fh:
        S = json.load(fh)
    runnable = sorted(S["runnable_T3plus"])
    tg = I.targets()
    tgd = {t["pdb"]: t for t in tg}
    nat = EV.native_torsions()
    inc = EV.incumbent()

    # ---- containment screen: no 6-mer of any target may occur in a corpus protein
    tset = set()
    for t in tg:
        s = t["seq"].upper()
        for i in range(len(s) - 5):
            tset.add(s[i:i + 6])
    leaks = 0
    seqs = {}
    for r in rows:
        seqs.setdefault(r["pdb"], []).append(r["aa"])
    for p, aa in seqs.items():
        s = "".join(aa)
        for i in range(len(s) - 5):
            if s[i:i + 6] in tset:
                leaks += 1
                break
    print("containment screen: {} / {} corpus proteins share a 6-mer with a target"
          .format(leaks, len(seqs)))

    ref, glob, sc = reference(rows)
    X, cphi, cpsi = corpus_xy(rows, ref, glob, sc)
    ycell = T.grid_bin(cphi, cpsi).astype(np.int64)
    groups = [r["pdb"] for r in rows]

    print("training the shift->torsion model (324-cell grid)")
    m_tor = fit(X, ycell, 324, "cls", groups=groups, tag="tors|shift")

    # sequence -> secondary-shift regressor, for the P-SHIFT arm
    print("training the sequence->shift model (for the P-SHIFT control)")
    Yn = np.full((len(rows), len(BACKBONE)), np.nan, np.float32)
    Xseq = X.copy()
    Xseq[:, :NF_SH] = 0.0                       # sequence context only
    for n, r in enumerate(rows):
        for k, a in enumerate(BACKBONE):
            v = r["sh"].get(a)
            if v is not None:
                Yn[n, k] = (v - ref.get((r["aa"], a), glob[a])) / sc[a]
    m_sh = fit(Xseq, Yn, len(BACKBONE), "reg", groups=groups, tag="shift|seq")

    # ------------------------------------------------------------------ evaluate arms
    seq_post = np.load(post_path("a_pepPos"))
    arms = ["E-SHIFT", "SHUF-SHIFT", "P-SHIFT", "REF-SHIFT", "MASK-ONLY", "SEQ-ONLY"]
    res = {a: {} for a in arms}
    errs = {a: [] for a in arms}
    posts = {a: {} for a in arms}

    for pid in runnable:
        t = tgd[pid]
        n = t["n"]
        aas = list(t["seq"])
        vals = target_shift_values(pid, A)
        gate = np.array(A[pid]["best"]["talos_mask"], bool)
        nphi, npsi = nat[pid]

        def build_feats(zero_values=False, force_present=False, use=None):
            src = vals if use is None else use
            return np.stack([feat_row(aas, src, i, ref, glob, sc,
                                      zero_values=zero_values,
                                      force_present=force_present)
                             for i in range(n)]).astype(np.float32)

        F = {}
        F["E-SHIFT"] = build_feats()
        F["REF-SHIFT"] = build_feats(zero_values=True, force_present=True)
        F["MASK-ONLY"] = build_feats(zero_values=True)
        # P-SHIFT: predict every secondary shift from sequence, then run the same model
        Xs = build_feats(zero_values=True, force_present=False)
        Xs[:, :NF_SH] = 0.0
        with torch.no_grad():
            pv = m_sh(torch.tensor(Xs)).numpy()
        pred_vals = []
        for i in range(n):
            d = {}
            for k, a in enumerate(BACKBONE):
                d[a] = float(pv[i, k]) * sc[a] + ref.get((aas[i], a), glob[a])
            pred_vals.append(d)
        F["P-SHIFT"] = build_feats(use=pred_vals)
        # SHUF-SHIFT: permute the measured shift dicts ACROSS residues of this target.  The
        # marginal shift distribution, the assignment pattern and the sequence are all
        # preserved; only the shift-to-residue correspondence -- the entire experimental
        # content -- is destroyed.  The strongest available null.
        rsh = np.random.default_rng(abs(hash(pid)) % (2 ** 31))
        F["SHUF-SHIFT"] = build_feats(use=[vals[j] for j in rsh.permutation(n)])

        for arm in arms:
            if arm == "SEQ-ONLY":
                P = seq_post[pid]
            else:
                P = posterior(m_tor, F[arm])
            posts[arm][pid] = P.astype(np.float32)
            d = T.decode(P)
            mp, ms = EV.determined_mask(n)
            res[arm].setdefault("sig", []).append(
                sigma_eff(d["phi"], d["psi"], nphi, npsi, (mp, ms)))
            errs[arm].append(T.wrap(d["psi"] - npsi))
            # direct build; residues failing the completeness gate are filled from the pool
            pf, pp, ca = EV.fill_pool(pid, d["phi"], d["psi"], gate)
            res[arm].setdefault("rmsd_gated", []).append(
                float(I.ca_rmsd(I.build_ca(pf, pp), ca)))
            res[arm].setdefault("rmsd_full", []).append(
                float(I.ca_rmsd(I.build_ca(d["phi"], d["psi"]), ca)))
            res[arm].setdefault("nll", []).append(
                float(np.mean(T.nll_native(P, nphi, npsi))))
            res[arm].setdefault("pdb", []).append(pid)

    print()
    print("=" * 92)
    print("REAL ARMS on the {} runnable targets.  incumbent on this subset = {:.4f} A"
          .format(len(runnable), float(np.mean([inc[p] for p in runnable]))))
    print("=" * 92)
    print("{:11s} {:>8s} {:>8s} {:>9s} {:>8s} {:>7s} {:>8s}".format(
        "arm", "sigma", "lag1", "emitted", "full", "<2A", "NLL"))
    incv = np.array([inc[p] for p in runnable])
    out = {"runnable": runnable, "incumbent_subset": float(incv.mean()),
           "incumbent_all126": float(np.mean(list(inc.values()))),
           "corpus_proteins": len(prots), "corpus_residues": len(rows),
           "containment_leaks": leaks}
    for arm in arms:
        sg = float(np.mean(res[arm]["sig"]))
        rg = np.array(res[arm]["rmsd_gated"])
        rf = np.array(res[arm]["rmsd_full"])
        l1 = lag1(errs[arm])
        print("{:11s} {:8.1f} {:8.3f} {:9.3f} {:8.3f} {:7.2f} {:8.3f}".format(
            arm, sg, l1, rg.mean(), rf.mean(), float((rg < 2.0).mean()),
            float(np.mean(res[arm]["nll"]))))
        pr = I.paired(rg, incv, folds=[tgd[p]["fold"] for p in runnable],
                      names=runnable)
        out[arm] = {"sigma_eff": sg, "lag1_error_autocorr": l1,
                    "emitted_gated": float(rg.mean()), "emitted_full": float(rf.mean()),
                    "frac_u2": float((rg < 2.0).mean()),
                    "nll": float(np.mean(res[arm]["nll"])),
                    "paired_vs_incumbent": pr,
                    "per_target": {p: float(v) for p, v in zip(runnable, rg)}}
    print()
    for arm in arms:
        pr = out[arm]["paired_vs_incumbent"]
        print("{:11s} paired vs incumbent: {}".format(arm, json.dumps(pr, default=str)[:220]))

    # --------- referencing diagnostic.  BMRB 13C/15N referencing offsets between labs are a
    # known hazard; a systematic per-entry offset would shift every feature.  Report the
    # per-target median SECONDARY shift by nucleus.  A well-referenced, conformationally
    # mixed entry sits near 0; a uniformly helical peptide legitimately does not, so this is
    # a diagnostic, not a correction (correcting it would remove real signal on a peptide
    # that is helical throughout).
    ref_diag = {}
    for nu in BACKBONE:
        med = []
        for pid in runnable:
            vals = target_shift_values(pid, A)
            aas = list(tgd[pid]["seq"])
            d = [vals[i][nu] - ref.get((aas[i], nu), glob[nu])
                 for i in range(len(aas)) if nu in vals[i]]
            if len(d) >= 4:
                med.append(float(np.median(d)))
        if med:
            m = np.array(med)
            ref_diag[nu] = {"n_targets": len(m), "median": float(np.median(m)),
                            "iqr": [float(np.percentile(m, 25)), float(np.percentile(m, 75))],
                            "frac_abs_gt_1ppm": float((np.abs(m) > 1.0).mean())}
    print()
    print("referencing diagnostic -- per-target median secondary shift (ppm):")
    for nu, v in ref_diag.items():
        print("  {:3s} n={:2d} median {:+.2f}  IQR [{:+.2f}, {:+.2f}]  |offset|>1ppm on {:.2f}"
              .format(nu, v["n_targets"], v["median"], v["iqr"][0], v["iqr"][1],
                      v["frac_abs_gt_1ppm"]))
    out["referencing_diagnostic"] = ref_diag

    cdir = os.path.join(ROOT, "s14", "cache")
    os.makedirs(cdir, exist_ok=True)
    for arm in arms:
        np.savez_compressed(os.path.join(cdir, "shift_post_%s.npz" % arm), **posts[arm])
    print("wrote per-arm 324-cell posteriors to s14/cache/shift_post_*.npz")

    out["what"] = ("REAL predictor arms. E-SHIFT uses MEASURED BMRB chemical shifts "
                   "(NMR-RESTRAINED, own column). P-SHIFT/REF-SHIFT/MASK-ONLY are nulls "
                   "with zero experimental content. Trained on an external BMRB+PDB corpus, "
                   "containment-screened against every target.")
    out["secs"] = round(time.time() - t0, 1)
    with open(os.path.join(RESULTS, "shift_model.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print("wrote s14/results/shift_model.json  [{:.0f}s]".format(time.time() - t0))


if __name__ == "__main__":
    main()
