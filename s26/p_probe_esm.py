"""s26/p_probe_esm.py -- lane P pre-sign-off PROBE. Reads NO native coordinate and NO RMSD.

Measures, under jobrun so the peak RSS is recorded honestly:
  (1) which compact ESM caches cover the distogram's training corpus and the 126 targets
      (esm_small.npz, s5/esm32.npz, s5/esmraw.npz, s7/repr_cache/esmcon.npz);
  (2) the per-fold training-set size (chains, pairs) of the shipped corpus
      (out-of-fold peptides + core.data.fold_fragments);
  (3) the resident size of s5/esmraw.npz once materialised as float32;
  (4) the wall time of ONE s12.instrument.project call on a PERSISTED production average
      (bench_results cache avg_ca; the comparison is against the persisted fit_ca / ca
      emissions, never the native).
"""
import os, sys, time, json
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)
import numpy as np, psutil
P = psutil.Process()
def rss(): return P.memory_info().rss / 1e9
out = {"rss_start_gb": rss()}
t0 = time.time()
from core import data as D
peps = D.load(); out["n_peptides"] = len(peps)
fc = getattr(D, "FRAGMENT_CACHE", None)
out["fragment_cache"] = fc; out["fragment_cache_exists"] = bool(fc and os.path.exists(fc))
frags = list(D.load_fragments()) if out["fragment_cache_exists"] else []
out["n_fragments"] = len(frags)
train_seqs = sorted({p.seq for p in peps} | {f.seq for f in frags})
out["n_train_seqs_unique"] = len(train_seqs)
from s12 import instrument as I
tg = I.targets(); tseqs = [t["seq"] for t in tg]; out["n_targets"] = len(tg)
folds = D.folds(5)
per_fold = {}
for f in range(5):
    ents = [p for p in peps if folds[p.seq] != f]
    fr = D.fold_fragments(f, 5) if frags else []
    npairs = sum(int(n * (n - 1) // 2 - (n - 1)) for n in (e.n for e in ents + fr))
    per_fold[f] = {"peptides": len(ents), "fragments": len(fr), "chains": len(ents) + len(fr), "pairs": npairs,
                   "targets_in_fold": int(sum(t["fold"] == f for t in tg))}
out["per_fold"] = per_fold; out["rss_after_db_gb"] = rss()
def cov(keys):
    ks = set(keys); c = sum(s in ks for s in train_seqs)
    return {"n_keys": len(ks), "train_covered": c, "train_missing": len(train_seqs) - c,
            "targets_covered": sum(s in ks for s in tseqs), "targets_missing": len(tseqs) - sum(s in ks for s in tseqs)}
z = np.load("esm_small.npz", allow_pickle=True); out["esm_small"] = cov([str(k) for k in z["seqs"]]); del z
out["rss_after_small_gb"] = rss()
z = np.load("s5/esm32.npz", allow_pickle=True); out["s5_esm32"] = cov([str(k) for k in z["keys"]]); del z
z = np.load("s7/repr_cache/esmcon.npz", allow_pickle=True); out["s7_esmcon"] = cov([str(k) for k in z["keys"]]); del z
t1 = time.time()
from s5 import esmraw
R = esmraw.load()
out["s5_esmraw"] = cov(R.keys()); out["s5_esmraw_load_s"] = round(time.time() - t1, 1); out["s5_esmraw_rss_gb"] = rss()
v = next(iter(R.values())); out["s5_esmraw_dim"] = list(v.shape); out["s5_esmraw_dtype"] = str(v.dtype)
out["s5_esmraw_total_MB_float32"] = round(sum(a.nbytes for a in R.values()) / 1e6, 1)
del R
rec = I.shipped_record(tg[0]["pdb"])
t2 = time.time()
pr = I.project(np.asarray(rec["avg_ca"], float), tg[0]["seq"], tg[0]["fold"])
out["project_target"] = tg[0]["pdb"]; out["project_wall_s"] = round(time.time() - t2, 1)
out["project_fit_ca_vs_persisted_maxabs"] = float(np.abs(pr["fit_ca"] - np.asarray(rec["fit_ca"], float)).max())
out["project_ca_vs_persisted_maxabs"] = float(np.abs(pr["ca"] - np.asarray(rec["ca"], float)).max())
out["peak_wset_gb"] = round(P.memory_info().peak_wset / 1e9, 3) if hasattr(P.memory_info(), "peak_wset") else None
out["wall_s"] = round(time.time() - t0, 1)
json.dump(out, open("s26/results/p_probe_esm.json", "w"), indent=1)
print(json.dumps(out, indent=1))
