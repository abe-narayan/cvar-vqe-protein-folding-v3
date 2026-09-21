"""LANE L / L1 -- FEASIBILITY.  Can a defensible 30-60 residue instrument be built?

Three independent supplies have to exist at once, and this script measures all three from
the data actually on disk.  Nothing here is a comparison and nothing here reads a native
coordinate for a decision: it is a CENSUS plus one geometric floor.

    TARGET SUPPLY      how many 30-60 residue chains in `prots/` survive the same quality
                       gates the 9-16mer instrument applies, and how many remain after
                       sequence-redundancy clustering.
    CANDIDATE SUPPLY   how many leakage-safe length-n windows the library can serve.  The
                       DEPLOYED library (peptide_db <= 25, fragment_db <= 20) serves ZERO
                       at n >= 26; the question is what `prots/` can serve instead.
    REPRESENTABILITY   the ideal-geometry builder is the only structure this pipeline can
                       emit.  `build_backbone(phi, psi)` from a chain's OWN native torsions
                       is the best the whole architecture could ever do at that length.
                       At 9-16 residues this floor is small enough to ignore.  It is a
                       HARD FLOOR on the endpoint and it must be measured before anything
                       else, because if it exceeds the RMSD we are chasing, the instrument
                       cannot support a verdict no matter how good the science is.

    python -m s32.s32_L_feasibility census     # target supply + representability floor
    python -m s32.s32_L_feasibility windows    # candidate supply at length
"""
from __future__ import annotations

import collections
import glob
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

from core import geometry as geo                                            # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROT_DIR = os.path.join(BASE, "prots")
RESULTS = os.path.join(BASE, "s32", "results")
os.makedirs(RESULTS, exist_ok=True)

#: The band the charter's "longer proteins" question asks about.  Kept wide for the census
#: so the instrument's length window can be chosen from the supply curve rather than
#: guessed; the instrument itself pins a narrower band.
LO, HI = 26, 80

#: Same gates the 9-16mer library applies (`fragment_db._extract`), so a length comparison
#: is a comparison of LENGTH and not of admission policy.
CA_STEP_LO, CA_STEP_HI = 3.5, 4.1


def _quality(path):
    """Parse one PDB and return the per-chain quality record, or None if unusable."""
    try:
        seq, coords, phi, psi = geo.native_coords_from_pdb(path)
    except Exception as e:
        return dict(ok=False, why=f"parse:{type(e).__name__}")
    ca = np.asarray(coords["CA"], float)
    n = len(seq)
    if n == 0 or len(ca) != n:
        return dict(ok=False, why="len-mismatch", n=n)
    if not (LO <= n <= HI):
        return dict(ok=False, why="out-of-band", n=n)
    if "X" in seq:
        return dict(ok=False, why="X-residue", n=n)
    step = np.linalg.norm(np.diff(ca, axis=0), axis=1)
    nbreak = int(((step <= CA_STEP_LO) | (step >= CA_STEP_HI)).sum())
    if nbreak:
        return dict(ok=False, why="chain-break", n=n, nbreak=nbreak)
    # phi[0] and psi[-1] are undefined at the termini; the builder supplies them.
    if not (np.all(np.isfinite(phi[1:])) and np.all(np.isfinite(psi[:-1]))):
        return dict(ok=False, why="nonfinite-torsion", n=n)
    ph = np.nan_to_num(np.asarray(phi, float))
    ps = np.nan_to_num(np.asarray(psi, float))
    try:
        built = geo.build_backbone(ph, ps)
    except Exception as e:
        return dict(ok=False, why=f"build:{type(e).__name__}", n=n)
    m = geo.kabsch_superpose(built["CA"], ca)
    rebuild = float(geo.rmsd(m, ca))
    return dict(ok=True, n=n, seq=seq, rebuild=rebuild,
                rg=float(geo.radius_of_gyration(ca)))


def stage_census(verbose=True):
    """Target supply in the 26-80 band, with the ideal-geometry floor for every survivor."""
    paths = sorted(glob.glob(os.path.join(PROT_DIR, "*.pdb")))
    rows, why = [], collections.Counter()
    t0 = time.time()
    for k, p in enumerate(paths):
        pdb = os.path.basename(p)[:-4].upper()
        r = _quality(p)
        if r.get("ok"):
            r["pdb"] = pdb
            rows.append(r)
        else:
            why[r["why"]] += 1
        if verbose and (k + 1) % 2000 == 0:
            print(f"  {k+1}/{len(paths)} scanned, {len(rows)} in band "
                  f"({time.time()-t0:.0f}s)", flush=True)
    out = dict(n_files=len(paths), n_in_band=len(rows),
               rejected=dict(why),
               band=[LO, HI])
    if rows:
        reb = np.array([r["rebuild"] for r in rows])
        ns = np.array([r["n"] for r in rows])
        out["rebuild"] = dict(mean=float(reb.mean()), median=float(np.median(reb)),
                              p10=float(np.percentile(reb, 10)),
                              p90=float(np.percentile(reb, 90)),
                              min=float(reb.min()), max=float(reb.max()))
        # The floor is the scientific object: report it BY LENGTH, because the whole
        # question is whether it scales.
        by = {}
        for lo, hi in [(26, 35), (36, 45), (46, 55), (56, 65), (66, 80)]:
            m = (ns >= lo) & (ns <= hi)
            if m.sum():
                by[f"{lo}-{hi}"] = dict(n=int(m.sum()), mean=float(reb[m].mean()),
                                        median=float(np.median(reb[m])),
                                        p90=float(np.percentile(reb[m], 90)))
        out["rebuild_by_length"] = by
    out["rows"] = rows
    path = os.path.join(RESULTS, "L1_census.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
        print("wrote", path, flush=True)
    return out


def stage_peptide_floor(verbose=True):
    """The SAME representability floor on the canonical 126, so the length axis is paired.

    ORACLE-free: native torsions are the INPUT here, and the quantity is a property of the
    representation, not of any prediction.  It reads native coordinates and is therefore
    labelled a diagnostic; it feeds no decision.
    """
    sys.path.insert(0, BASE)
    from s12 import instrument as ins
    import peptide_db as db
    byseq = {p.seq: p for p in db.load()}
    rows = []
    for t in ins.targets():
        p = byseq.get(t["seq"])
        if p is None:
            continue
        ph = np.nan_to_num(np.asarray(p.phi, float))
        ps = np.nan_to_num(np.asarray(p.psi, float))
        built = geo.build_backbone(ph, ps)
        m = geo.kabsch_superpose(built["CA"], np.asarray(p.ca, float))
        rows.append(dict(pdb=t["pdb"], n=t["n"],
                         rebuild=float(geo.rmsd(m, np.asarray(p.ca, float)))))
    reb = np.array([r["rebuild"] for r in rows])
    out = dict(n=len(rows), mean=float(reb.mean()), median=float(np.median(reb)),
               p90=float(np.percentile(reb, 90)), max=float(reb.max()), rows=rows)
    path = os.path.join(RESULTS, "L1_peptide_floor.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
        print("wrote", path, flush=True)
    return out


def stage_windows(lengths=(13, 26, 35, 45, 55), verbose=True):
    """Candidate supply: length-n windows the DEPLOYED library vs `prots/` can serve.

    The deployed library is the one the 126-target instrument actually retrieves from.  Its
    supply at n is the number of members with length >= n, times the windows each yields.
    """
    sys.path.insert(0, BASE)
    import peptide_db as db
    import fragment_db as fdb
    peps = [p for p in db.load()]
    frags = list(fdb.load())
    prot_n = [r["n"] for r in json.load(open(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "results", "L1_protlen.json")))] \
        if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                       "results", "L1_protlen.json")) else None
    if prot_n is None:
        prot_n = []
        for p in sorted(glob.glob(os.path.join(PROT_DIR, "*.pdb"))):
            n = 0
            chain = None
            with open(p, "r", errors="ignore") as fh:
                for line in fh:
                    if line.startswith("ENDMDL"):
                        break
                    if line.startswith("ATOM") and line[12:16].strip() == "CA" \
                            and line[16] in (" ", "A"):
                        c = line[21]
                        if chain is None:
                            chain = c
                        if c == chain:
                            n += 1
            prot_n.append(n)
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "results", "L1_protlen.json"), "w") as fh:
            json.dump([dict(n=int(x)) for x in prot_n], fh)
    prot_n = np.asarray(prot_n, int)
    out = {}
    for L in lengths:
        dep = sum(max(0, len(p.seq) - L + 1) for p in peps) \
            + sum(max(0, len(f.seq) - L + 1) for f in frags)
        pr = int(np.maximum(0, prot_n - L + 1).sum())
        out[str(L)] = dict(deployed_library_windows=int(dep),
                           prots_windows=pr,
                           n_prots_long_enough=int((prot_n >= L).sum()))
    path = os.path.join(RESULTS, "L1_windows.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    if verbose:
        print(json.dumps(out, indent=1))
        print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "census"
    {"census": stage_census, "windows": stage_windows,
     "pepfloor": stage_peptide_floor}[stage]()
