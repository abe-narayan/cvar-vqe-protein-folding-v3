# QUARANTINE — PROOF-BUILD OUTPUT. NOT RESULTS. DO NOT PUBLISH, CITE, OR RESTORE.

Quarantined 2026-09-08 by the Sprint 25 coordinator on the audit lane's report.

## What this is

1,134 PDB files, 126 under each of the seven real comparison configuration names plus
`production` and `native`, produced by `--mode proof`. **Six of the seven comparison arms are
`native + N(0, sigma)` with sigma 1.3-2.2.** They are not predictions.

## Why it was dangerous

- **Not one of the 1,134 files carried any synthetic or provenance marker.** `is_synthetic`
  propagated into the JSON payload but never reached the PDB headers, whose ten REMARK 999 keys
  are BASIS, CONFIGURATION, FOLD, LENGTH, PDB_ID, REFERENT, RMSD_CA, RMSD_CA_MEMORY, SEQUENCE,
  TARGET_ID. **Separated from `results.json`, each file is indistinguishable from a real result.**
- `leaderboard.json` ranked `legacy_amber_distogram` first at **2.0921 A**, with a paired effect of
  **-0.956 A against production at 3.0483** — a fake method beating the real one, in a file that
  reads as a finished deliverable.
- **4 of the 126 `distogram` files are EXACT NATIVES, `RMSD_CA 0.000000`** — test fixtures, because
  `test_export.py` wrote into the real output directory under real configuration names.
- `T006__amber.pdb` carries `RMSD_CA PENDING`: `export_prediction` writes twice and is not atomic.
- The adopted pool-oracle guard (refuse mean RMSD < pool_best = 1.7108 A) **would not have caught any
  of it** — every synthetic arm's mean is above the ceiling by construction.

## What was NOT touched

`results/benchmark_manifest.json` and `results/monomer_manifest.json` are git-tracked project files
predating this sprint. `results/summary/target_map.json` is a pure T-number-to-PDB mapping and
carries no RMSD; both remain in place.

Retained as evidence for the final report's leakage-controls section. Regenerate nothing from here.
