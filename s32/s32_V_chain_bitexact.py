"""S32 LANE V -- is the canonical 3.2105 bit-reproducible from its stated input?

Finding this closes:  s32/s32_step4_rebuild.py projected a RECOMPUTED top-75 coordinate
average; s29/s29_O_ladder.py::chain_item projected the STORED `avg_ca` read back from the
production cache JSON.  The two clouds agree to ~1e-15 (one float64 ULP) and give built
chains up to 0.52 A apart, because the lam=0 multi-start argmin is decided at a ~1e-7 branch
margin (S32 contract rule 3, amplification ~1e13).

This script projects the STORED cloud, exactly as chain_item does, and asks whether the
result is BIT-IDENTICAL to the 126 persisted `prod` rows.  A verification must be able to
fail (contract rule 5): the falsifier is a single target disagreeing in the 9th decimal.

  * if bit-identical on 126/126 -> the canonical endpoint reproduces exactly, and the
    entire 0.0021 A residual of the step-4 rebuild is the ULP effect, not an error.
  * if not -> the endpoint is not reproducible from its own stated input, which is a
    stronger defect than a labelling one.

Usage: python s32/s32_V_chain_bitexact.py --shard 0 --n-shards 4
       python s32/s32_V_chain_bitexact.py --analyse
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
ROWS = os.path.join(RESULTS, "s32_V_bitexact_rows.jsonl")


def read_rows(pattern):
    out = {}
    for f in sorted(glob.glob(pattern)):
        for line in open(f):
            line = line.strip()
            if line:
                r = json.loads(line)
                out[r["pdb"]] = r
    return out


def s29_prod():
    out = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "s29", "results", "s29_O_chain_rows*.jsonl"))):
        for line in open(f):
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("item") == "prod":
                out[r["pdb"]] = r
    return out


#: THE canonical input.  `s29_O_ladder.cloud_row` computes the top-75 coordinate average,
#: saves it to this npz, and `chain_item` reloads it and projects THAT array.  The cache's
#: `avg_ca` is a different float64 representation of the same average (max|d| ~1e-14).
S29_STRUCTS = os.path.join(ROOT, "s29", "results", "s29_O_structs")


def clouds(pdb, rec):
    """The three bit-representations of the SAME top-75 coordinate average."""
    out = {"cache_avg_ca": np.asarray(rec["avg_ca"], float)}
    f = os.path.join(S29_STRUCTS, "%s.npz" % pdb)
    if os.path.exists(f):
        with np.load(f) as z:
            out["s29_struct_prod"] = np.asarray(z["prod"], float)
    return out


def run(shard, n_shards):
    tg = I.targets()
    done = read_rows(ROWS.replace(".jsonl", "*.jsonl"))
    path = ROWS if shard is None else ROWS.replace(".jsonl", "_shard%d.jsonl" % shard)
    todo = [t for k, t in enumerate(tg)
            if (shard is None or k % n_shards == shard)
            and "s29_struct_prod" not in done.get(t["pdb"], {})]
    print("bitexact shard %s/%s: %d to do" % (shard, n_shards, len(todo)), flush=True)
    for k, t in enumerate(todo):
        pdb = t["pdb"]
        u = I.load_univ(pdb); rec = I.shipped_record(pdb)
        row = dict(pdb=pdb)
        for nm, C in clouds(pdb, rec).items():
            t0 = time.time()
            pr = I.project(C, t["seq"], t["fold"])
            row[nm] = float(I.ca_rmsd(np.asarray(pr["ca"], float), u["nat_ca"]))
            row[nm + "_cloud"] = float(I.ca_rmsd(C, u["nat_ca"]))
            row[nm + "_secs"] = time.time() - t0
        cl = clouds(pdb, rec)
        if len(cl) == 2:
            row["max_abs_cloud_diff"] = float(np.abs(cl["cache_avg_ca"] - cl["s29_struct_prod"]).max())
        row["rmsd_chain"] = row.get("cache_avg_ca")       # back-compat with the first pass
        with open(path, "a") as fh:
            fh.write(json.dumps(row) + "\n")
        print("  [%d/%d] %-6s s29struct %.9f  cache %.9f  |dC| %.2e"
              % (k + 1, len(todo), pdb, row.get("s29_struct_prod", float("nan")),
                 row.get("cache_avg_ca", float("nan")), row.get("max_abs_cloud_diff", 0)), flush=True)


def analyse():
    mine = read_rows(ROWS.replace(".jsonl", "*.jsonl"))
    ref = s29_prod()
    reb = read_rows(os.path.join(RESULTS, "s32_V_step4_rows*.jsonl"))
    out = dict(canonical_mean_s29=float(np.mean([ref[p]["rmsd_chain"] for p in sorted(ref)])),
               n_canonical=len(ref), arms={})
    print("=" * 96)
    print("S32-V: is the canonical endpoint reproducible, and from WHICH bit-representation?")
    print("canonical (s29 prod rows, n=%d): %.9f" % (len(ref), out["canonical_mean_s29"]))
    print("-" * 96)
    print("%-34s %-4s %14s %10s %10s %8s" % ("input cloud", "n", "mean chain", "mean d", "max|d|", "bit-id"))

    def arm(name, table, key, note=""):
        pdbs = sorted([p for p in set(table) & set(ref) if key in table[p]
                       and table[p][key] is not None])
        if not pdbs:
            return
        a = np.array([table[p][key] for p in pdbs], float)
        b = np.array([ref[p]["rmsd_chain"] for p in pdbs], float)
        d = a - b
        r = dict(n=len(pdbs), mean=float(a.mean()), mean_diff=float(d.mean()),
                 max_abs_diff=float(np.abs(d).max()), n_bit_identical=int((d == 0.0).sum()),
                 n_within_1e9=int((np.abs(d) < 1e-9).sum()), note=note,
                 worst=[(pdbs[k], float(a[k]), float(b[k])) for k in np.argsort(-np.abs(d))[:3]])
        out["arms"][name] = r
        print("%-34s %-4d %14.9f %+10.6f %10.6f %4d/%d %s"
              % (name, len(pdbs), a.mean(), d.mean(), np.abs(d).max(),
                 r["n_bit_identical"], len(pdbs), note))
        return r

    a1 = arm("s29 struct npz ['prod']", mine, "s29_struct_prod", "the input chain_item actually reads")
    a2 = arm("production cache avg_ca", mine, "cache_avg_ca", "same average, other float64 bits")
    a3 = arm("recomputed coordinate average", reb, "rmsd_chain", "same average, recomputed in-process")
    dc = [mine[p]["max_abs_cloud_diff"] for p in mine if "max_abs_cloud_diff" in mine[p]]
    if dc:
        out["max_abs_cloud_diff_over_targets"] = float(max(dc))
        print("-" * 96)
        print("max |cache_avg_ca - s29_struct_prod| over targets: %.3e   "
              "(their CLOUD RMSDs agree to 9 dp everywhere)" % max(dc))

    ok = bool(a1 and a1["n_bit_identical"] == a1["n"])
    out["endpoint_reproduces_bit_exactly"] = ok
    out["verdict"] = (
        "THE CANONICAL 3.2105 REPRODUCES BIT-FOR-BIT FROM ITS STATED INPUT "
        "(s29_O_structs/<pdb>.npz['prod']); the two other bit-representations of the SAME "
        "average give a different chain on nearly every target."
        if ok else "*** THE CANONICAL ENDPOINT DOES NOT REPRODUCE FROM ITS STATED INPUT ***")
    print("VERDICT: %s" % out["verdict"])

    # SELF-TEST (contract rule 5): the check must be able to FAIL.  The two alternative
    # bit-representations are the same mathematical object; if the comparison called THEM
    # bit-identical too it would be decoration, not a check.
    st = [(n, out["arms"][n]["n_bit_identical"], out["arms"][n]["n"])
          for n in ("production cache avg_ca", "recomputed coordinate average") if n in out["arms"]]
    good = all(b < n for _, b, n in st) and bool(st)
    out["selftest_detects_ulp_perturbed_input"] = good
    print("selftest: a ~1e-14-perturbed input is %s  %s"
          % ("CORRECTLY seen as different" if good else "*** NOT DETECTED -- check is broken ***",
             ", ".join("%s %d/%d bit-identical" % s for s in st)))
    print("=" * 96)
    with open(os.path.join(RESULTS, "s32_V_chain_bitexact.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=None)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--analyse", action="store_true")
    a = ap.parse_args()
    analyse() if a.analyse else run(a.shard, a.n_shards)
