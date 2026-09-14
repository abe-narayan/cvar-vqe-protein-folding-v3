# verify/ audits re-run under the S26 governor (lane I)

Each audit was executed in-process by `s26/i_verify_rerun.py run <name>` with its writes redirected to `s26/results/verify/`; the tracked JSON's sha256 was asserted unchanged. Leaves compared: identical / close (float noise) / different / run-property (clocks, hosts, commits) / added / removed.

| audit | tag | peak RSS | wall | tracked JSON | verdict | identical | close | different | run-prop | added | removed | error |
|---|---|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---|
| vqe_lfo_audit | CPU | 0.002 GB | 0.3 s | vqe_lfo_audit.json | IDENTICAL | 32 | 0 | 0 | 0 | 0 | 0 |  |
| cvar_audit | CPU | 0.003 GB | 0.9 s | cvar_audit.json | IDENTICAL | 18 | 0 | 0 | 0 | 0 | 0 |  |
| ansatz_audit | CPU | 0.005 GB | 3.3 s | ansatz_audit.json | IDENTICAL | 19 | 0 | 0 | 0 | 0 | 0 |  |
| headline_audit | CPU | 0.002 GB | 0.6 s | headline_audit.json | DIFFERS | 123 | 0 | 0 | 0 | 8 | 0 |  |
| legacy_audit | CPU | 0.127 GB | 27.8 s | legacy_audit.json | IDENTICAL | 30 | 0 | 0 | 0 | 0 | 0 |  |
| amber_audit | AMBER | 0.265 GB | 42.2 s | amber_audit.json | IDENTICAL | 31 | 0 | 0 | 0 | 0 | 0 |  |
| amber_platform | AMBER | 0.304 GB | 70.4 s | amber_platform.json | NO FRESH OUTPUT |  |  |  |  |  |  |  |
| leak_audit | AMBER | 0.581 GB | 118.0 s | leak_audit.json | DIFFERS | 24 | 0 | 1 | 0 | 0 | 0 |  |
| grad_key_collision | AMBER | 0.578 GB | 2995.1 s | grad_key_collision.json | DIFFERS | 7 | 0 | 3 | 26 | 0 | 10 |  |
| equiv_compare | CPU | 0.003 GB | 0.8 s | equivalence.json | DIFFERS | 924 | 0 | 5 | 0 | 1 | 0 |  |
| projection_divergence | CPU | 0.002 GB | 1.2 s | projection_divergence.json | NO FRESH OUTPUT |  |  |  |  |  |  | ValueError: shape mismatch: () vs (12, 3) |
| project_arms | CPU | 0.008 GB | 0.3 s | project_arms.json | DIFFERS | 47 | 0 | 124 | 0 | 0 | 0 |  |

## Not re-runnable

- `verify/recon_containment_audit.json`: written by a deleted script (leak_audit_mine.json, moved in the 2026-09-04 consolidation); no runner in verify/
- `verify/project_inputs.json`: written by `core.project harvest` from the production cache; it is the INPUT the projection audits read, re-harvesting it is covered by project_exactness / project_equiv reading it

## First differences per audit

### leak_audit / leak_audit.json
- `backends.project`: tracked `"s8.project"`, fresh `"core.project"`

### grad_key_collision / grad_key_collision.json
- `COLLISION`: tracked `true`, fresh `false`
- `config_has_a_field_for_the_gradient`: tracked `false`, fresh `true`
- `science_fields`: tracked `["amber", "amber_k", "amber_steps", "dev_mode", "k", "lam", "legacy", "legacy_top", "m", "maxiter", "min_sep", "multi_start", "n_folds", "penalty", "quantum"...`, fresh `["amber", "amber_k", "amber_steps", "dev_mode", "k", "lam", "legacy", "legacy_top", "m", "maxiter", "min_sep", "multi_start", "n_folds", "penalty", "project_...`

### equiv_compare / equivalence.json
- `baseline_dir`: tracked `"C:\\Users\\abena\\Protein-Folding-Algorithm\\bench_results\\cache\\29cc4e7cba2feeff"`, fresh `"C:\\Users\\abena\\Protein-Folding-Algorithm\\bench_results\\cache\\44a9305c80808776"`
- `baseline_key`: tracked `"29cc4e7cba2feeff"`, fresh `"44a9305c80808776"`
- `n_optimised`: tracked `8`, fresh `126`
- `optimised_dir`: tracked `"C:\\Users\\abena\\Protein-Folding-Algorithm\\bench_results\\cache\\40774ae7b037d259"`, fresh `"C:\\Users\\abena\\Protein-Folding-Algorithm\\bench_results\\cache\\1fc9f2dcf489e2fb"`
- `optimised_key`: tracked `"40774ae7b037d259"`, fresh `"1fc9f2dcf489e2fb"`

### project_arms / project_arms.json
- `consolidation_faithful__baseline_vs_opt_fd.amber_ca.bit_identical`: tracked `"0/8"`, fresh `"119/119"`
- `consolidation_faithful__baseline_vs_opt_fd.amber_ca.worst_max_abs`: tracked `22.333503366363715`, fresh `0.0`
- `consolidation_faithful__baseline_vs_opt_fd.avg_ca.bit_identical`: tracked `"8/8"`, fresh `"126/126"`
- `consolidation_faithful__baseline_vs_opt_fd.ca.bit_identical`: tracked `"0/8"`, fresh `"126/126"`
- `consolidation_faithful__baseline_vs_opt_fd.ca.worst_max_abs`: tracked `22.292044613202613`, fresh `0.0`
- `consolidation_faithful__baseline_vs_opt_fd.fit_ca.bit_identical`: tracked `"0/8"`, fresh `"126/126"`
- `consolidation_faithful__baseline_vs_opt_fd.fit_ca.worst_max_abs`: tracked `27.111554683975196`, fresh `0.0`
- `consolidation_faithful__baseline_vs_opt_fd.n`: tracked `8`, fresh `126`
- `consolidation_faithful__baseline_vs_opt_fd.phi.bit_identical`: tracked `"0/8"`, fresh `"126/126"`
- `consolidation_faithful__baseline_vs_opt_fd.phi.worst_max_abs`: tracked `6.285326248455182`, fresh `0.0`
- `consolidation_faithful__baseline_vs_opt_fd.pool_best.bit_identical`: tracked `"8/8"`, fresh `"126/126"`
- `consolidation_faithful__baseline_vs_opt_fd.pool_best.mean_a`: tracked `1.4937286227941513`, fresh `1.7108244199364904`

