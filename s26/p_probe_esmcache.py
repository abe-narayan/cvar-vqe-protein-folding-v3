"""s26/p_probe_esmcache.py -- the ONE permitted load of esm_cache.npz (1.5 GB object npz), under
jobrun --tag ESM, to record what it costs in RAM and time. Loads, counts, exits. No structure read."""
import os, sys, time, json
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
import numpy as np, psutil
P = psutil.Process(); t0 = time.time()
z = np.load("esm_cache.npz", allow_pickle=True)
seqs = z["seqs"]; t_seqs = time.time() - t0; r_seqs = P.memory_info().rss / 1e9
vals = z["vals"]; t_vals = time.time() - t0; r_vals = P.memory_info().rss / 1e9
n = int(len(seqs)); step = max(1, n // 500); dims = set(); nres = 0
for k in range(0, n, step):
    a = np.asarray(vals[k][0]); dims.add(int(a.shape[1])); nres += a.shape[0]
out = {"n_seqs": n, "t_seqs_s": round(t_seqs, 1), "rss_after_seqs_gb": round(r_seqs, 3),
       "t_vals_s": round(t_vals, 1), "rss_after_vals_gb": round(r_vals, 3),
       "rep_dims_sampled": sorted(dims), "mean_len_sampled": round(nres / max(1, len(range(0, n, step))), 2),
       "peak_wset_gb": round(P.memory_info().peak_wset / 1e9, 3) if hasattr(P.memory_info(), "peak_wset") else None,
       "wall_s": round(time.time() - t0, 1)}
json.dump(out, open("s26/results/p_probe_esmcache.json", "w"), indent=1)
print(json.dumps(out, indent=1))
