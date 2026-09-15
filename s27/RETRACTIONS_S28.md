# RETRACTIONS_S28 (append-only; lane D keeps it)

Every S28 retraction, with the ledger entry it retracts, the entry that retracts it, and the
artefact. Nothing superseded is deleted anywhere; this file is the index.

| id | retracted claim (entry) | retracting entry | artefact | disposition |
|---|---|---|---|---|
| R1 | S28-L7 (lane C): "built-chain contrasts carry a numerical floor near 1e-5 A" (one target) | S28-L18 (lane D) | `s27/results/s28_A_chain_rows.jsonl` vs `s27/results/chain_rows.jsonl :: DIS`, 22 targets: max 0.0186, mean 0.0025 | SCOPE CORRECTION: the floor is 1e-5 to 2e-2 per target (0.003 mean) whenever the two sides of a chain contrast come from clouds that differ at floating-point level; zero when they share a code path. No result changes. |
| R2 | S28-L18 (lane D): the built chain's floor is "0.02 A per target (0.0025 mean)" (22 targets) | S28-L27 (lane D) | `s27/results/s28_A_chain_rows.jsonl :: prod` vs `chain_rows.jsonl :: DIS`, 126 targets: mean -0.0056, 12 above 0.02, 3 above 0.1, 2LNG -0.513 | SCOPE CORRECTION (my own): the tail reaches 0.5 A on one target in 126 and the mean floor is 0.006 A; the reading (branch flips under floating-point input differences; two sides of a chain contrast must share a code path) stands. No result changes. |
