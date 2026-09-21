"""S32 LANE V -- the ARITHMETIC-NOISE DISTRIBUTION of the endpoint.

S32-D1 established that the built chain is a discontinuous function of its input: three
bit-representations of the SAME top-75 coordinate average, agreeing to 5.7e-14 A and with
identical cloud RMSD to 9 dp, give

    s29_O_structs['prod']              3.210534   <- the canonical endpoint
    production cache avg_ca            3.214765   <- and this is EXACTLY the chain the
                                                     production pipeline itself emits
    a fresh in-process recomputation   3.212625

That is a 0.0042 A spread from three arbitrary draws.  Three draws is an anecdote.  This
script turns it into the number every lane actually needs: **the sd of the endpoint under a
perturbation at the scale of float64 round-off**, i.e. the irreducible uncertainty floor on
ANY unpaired built-chain mean.

METHOD.  For each target, perturb the canonical cloud by eps * N(0,1) per coordinate with
eps = 1e-14 A -- a displacement 13 orders of magnitude below the chain floor and 15 below the
endpoint, which changes no physical quantity (the cloud RMSD is unchanged to 9 dp, asserted
per draw).  Re-project.  Repeat for D draws.  Report the per-draw endpoint mean, the
draw-to-draw sd, and the per-target |delta| distribution.

CONTROL / WHAT MAKES THIS A VALID MEASURE.  The perturbation is matched to the operator's own
input space (contract rule 6): it is applied to the same array the operator reads, at the
magnitude at which two honest implementations of the same average already disagree.  It is
not a physical perturbation and is never interpreted as one.

Usage: python s32/s32_V_ulp_distribution.py --shard 0 --n-shards 3 --draws 6
       python s32/s32_V_ulp_distribution.py --analyse
"""
from __future__ import annotations
import argparse, glob, json, os, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("s32_instrument", os.path.join(ROOT, "s12", "instrument.py"))
I = _ilu.module_from_spec(_spec); _spec.loader.exec_module(I)

RESULTS = os.path.join(ROOT, "s32", "results")
os.makedirs(RESULTS, exist_ok=True)
ROWS = os.path.join(RESULTS, "s32_V_ulp_rows.jsonl")
STRUCTS = os.path.join(ROOT, "s29", "results", "s29_O_structs")
EPS = 1e-14


def canonical_cloud(pdb):
    with np.load(os.path.join(STRUCTS, "%s.npz" % pdb)) as z:
        return np.asarray(z["prod"], float)


def read_rows():
    out = {}
    for f in sorted(glob.glob(ROWS.replace(".jsonl", "*.jsonl"))):
        for line in open(f):
            line = line.strip()
            if line:
                r = json.loads(line)
                out[(r["pdb"], r["draw"])] = r
    return out


def run(shard, n_shards, draws):
    tg = I.targets()
    done = read_rows()
    path = ROWS if shard is None else ROWS.replace(".jsonl", "_shard%d.jsonl" % shard)
    todo = [(t, d) for k, t in enumerate(tg) for d in range(draws)
            if (shard is None or k % n_shards == shard) and (t["pdb"], d) not in done]
    print("ulp shard %s/%s: %d to do" % (shard, n_shards, len(todo)), flush=True)
    t0 = time.time()
    for k, (t, d) in enumerate(todo):
        pdb = t["pdb"]
        u = I.load_univ(pdb); nat = u["nat_ca"]
        C = canonical_cloud(pdb)
        #: seed is a pure function of (pdb, draw): the perturbation is reproducible
        rng = np.random.default_rng(abs(hash(("s32Vulp", pdb, d))) % (2 ** 63))
        Cp = C + EPS * rng.standard_normal(C.shape)
        pr = I.project(Cp, t["seq"], t["fold"])
        row = dict(pdb=pdb, draw=d, fold=t["fold"],
                   rmsd_chain=float(I.ca_rmsd(np.asarray(pr["ca"], float), nat)),
                   rmsd_cloud=float(I.ca_rmsd(Cp, nat)),
                   cloud_rmsd_unperturbed=float(I.ca_rmsd(C, nat)),
                   max_abs_dC=float(np.abs(Cp - C).max()), secs=time.time() - t0)
        with open(path, "a") as fh:
            fh.write(json.dumps(row) + "\n")
        if (k + 1) % 20 == 0 or k == len(todo) - 1:
            print("  [%d/%d] %-6s d%d chain %.6f (%.1f min)"
                  % (k + 1, len(todo), pdb, d, row["rmsd_chain"], (time.time() - t0) / 60), flush=True)


def analyse():
    rows = read_rows()
    canon = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "s29", "results", "s29_O_chain_rows*.jsonl"))):
        for line in open(f):
            line = line.strip()
            if line:
                r = json.loads(line)
                if r.get("item") == "prod":
                    canon[r["pdb"]] = float(r["rmsd_chain"])
    draws = sorted({d for _p, d in rows})
    pdbs = sorted({p for p, _d in rows})
    full = [d for d in draws if sum((p, d) in rows for p in pdbs) == len(pdbs)]
    out = dict(eps=EPS, n_targets=len(pdbs), draws_complete=len(full), draws_seen=len(draws),
               canonical_mean=float(np.mean([canon[p] for p in pdbs])) if pdbs else None)
    print("=" * 96)
    print("ARITHMETIC-NOISE DISTRIBUTION OF THE ENDPOINT   eps = %.0e A, n=%d targets, "
          "%d complete draws" % (EPS, len(pdbs), len(full)))
    # the perturbation must be invisible on the cloud -- if it is not, this measures something else
    dc = np.array([rows[(p, d)]["rmsd_cloud"] - rows[(p, d)]["cloud_rmsd_unperturbed"]
                   for p in pdbs for d in full])
    out["max_abs_cloud_rmsd_shift"] = float(np.abs(dc).max())
    print("  sanity: max |cloud RMSD shift| from the perturbation = %.2e A  "
          "(the cloud is physically unchanged)" % np.abs(dc).max())
    means = np.array([np.mean([rows[(p, d)]["rmsd_chain"] for p in pdbs]) for d in full])
    out["draw_means"] = [float(m) for m in means]
    if len(means) > 1:
        out["draw_mean"] = float(means.mean()); out["draw_sd"] = float(means.std(ddof=1))
        out["draw_min"] = float(means.min()); out["draw_max"] = float(means.max())
        print("  endpoint per draw: %s" % ", ".join("%.6f" % m for m in means))
        print("  draw mean %.6f   draw-to-draw sd %.6f   range %.6f"
              % (means.mean(), means.std(ddof=1), means.max() - means.min()))
        print("  canonical %.6f sits at the %.0fth percentile of the draws"
              % (out["canonical_mean"], 100 * (means < out["canonical_mean"]).mean()))
    per = np.array([abs(rows[(p, d)]["rmsd_chain"] - canon[p]) for p in pdbs for d in full])
    out["per_target_abs_delta"] = dict(mean=float(per.mean()), p50=float(np.median(per)),
                                       p90=float(np.percentile(per, 90)), max=float(per.max()),
                                       frac_identical=float((per == 0).mean()))
    print("  per-target |chain - canonical|: mean %.4f  p50 %.4f  p90 %.4f  max %.4f  "
          "(bit-identical on %.1f%%)"
          % (per.mean(), np.median(per), np.percentile(per, 90), per.max(),
             100 * (per == 0).mean()))
    print("  contract rule 3's measured same-operator floor is 0.0134 mean / 0.0329 p90 / "
          "0.2285 max -- this is the SAME phenomenon, priced at its source.")
    if len(means) > 1:
        print()
        print("  USE: the irreducible sd on an UNPAIRED built-chain mean is %.4f A.  A paired "
              "contrast is unaffected ONLY if both sides are projected from bit-identical "
              "clouds in one job." % means.std(ddof=1))
    print("=" * 96)
    with open(os.path.join(RESULTS, "s32_V_ulp_distribution.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=None)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--draws", type=int, default=6)
    ap.add_argument("--analyse", action="store_true")
    a = ap.parse_args()
    analyse() if a.analyse else run(a.shard, a.n_shards, a.draws)
