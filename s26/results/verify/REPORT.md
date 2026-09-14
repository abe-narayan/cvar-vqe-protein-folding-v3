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
| determinism_audit | AMBER | 0.579 GB | 59.8 s | _ | no tracked JSON to diff |  |  |  |  |  |  |  |
| hazard_audit | AMBER | 0.525 GB | 21.7 s | _ | no tracked JSON to diff |  |  |  |  |  |  |  |
| project_selfcheck | CPU | 0.003 GB | 0.6 s | _ | no tracked JSON to diff |  |  |  |  |  |  |  |
| project_exactness | CPU | 0.092 GB | 587.0 s | project_exactness.json | DIFFERS | 1894 | 0 | 1 | 1 | 2 | 0 |  |
| project_equiv | CPU | 0.275 GB | 1798.9 s | project_equiv.json | DIFFERS | 12431 | 0 | 508 | 10 | 2 | 0 |  |
| project_degeneracy | CPU | 0.09 GB | 239.5 s | project_degeneracy.json | DIFFERS | 256 | 1 | 510 | 0 | 2 | 0 |  |
| project_iters | CPU | 0.273 GB | 384.7 s | project_iters.json | DIFFERS | 4 | 0 | 14 | 4 | 7 | 0 |  |
| project_stability | CPU | 0.273 GB | 1308.6 s | project_stability_partial77.json | DIFFERS | 469 | 0 | 2 | 0 | 0 | 1 |  |

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

### project_exactness / project_exactness.json
- `agg.s_per_target`: tracked `4.568959217458975`, fresh `4.650387336503223`

### project_equiv / project_equiv.json
- `agg.s_per_target.ref`: tracked `10.233509180953703`, fresh `7.689821138097498`
- `agg.s_per_target.ex`: tracked `4.674835245236792`, fresh `3.8137957515889807`
- `agg.s_per_target.fd`: tracked `2.367491449206321`, fresh `2.1160124063506602`
- `agg.s_per_target.an`: tracked `0.6235212849250751`, fresh `0.61906772143104`
- `rows[0].t_ref`: tracked `7.243612399965059`, fresh `7.445828400028404`
- `rows[0].t_ex`: tracked `3.496815799968317`, fresh `3.9513764000148512`
- `rows[0].t_fd`: tracked `2.352735499967821`, fresh `2.6301505000446923`
- `rows[0].t_an`: tracked `0.5057675000280142`, fresh `0.9298726000124589`
- `rows[1].t_ref`: tracked `9.680596000049263`, fresh `9.970852200000081`
- `rows[1].t_ex`: tracked `4.774421499983873`, fresh `5.054604000004474`
- `rows[1].t_fd`: tracked `2.4140376000432298`, fresh `3.0481882000458427`
- `rows[1].t_an`: tracked `0.5553847000119276`, fresh `0.9945384999737144`

### project_degeneracy / project_degeneracy.json
- `agg.gap_median`: tracked `0.06665769745825689`, fresh `0.06428815619227479`
- `agg.gap_min`: tracked `1.1851245765998053e-09`, fresh `2.804459381736124e-09`
- `agg.n_gap_under_1e4`: tracked `23`, fresh `22`
- `agg.runner_up_dist_median`: tracked `0.84600201106656`, fresh `0.8767447079132815`
- `agg.runner_up_dist_max`: tracked `2.6061174526946895`, fresh `2.468616994488777`
- `agg.n_tied_1e3_far_0p5`: tracked `1`, fresh `3`
- `agg.n_tied_1e2_far_0p5`: tracked `8`, fresh `10`
- `rows[0].best`: tracked `0.5195771787229689`, fresh `0.5197169327647683`
- `rows[0].gap_to_runner_up`: tracked `0.0011067270120530548`, fresh `0.04292251081244769`
- `rows[0].spread`: tracked `1.8776938351082784`, fresh `1.876491384064896`
- `rows[0].runner_up_dist`: tracked `0.03783408552742756`, fresh `0.6583542333408575`
- `rows[1].best`: tracked `1.2481900457886193`, fresh `1.2478925901618894`

### project_iters / project_iters.json
- `arms[0].solves`: tracked `90`, fresh `216`
- `arms[0].nit_mean`: tracked `161.33333333333334`, fresh `181.06944444444446`
- `arms[0].nfev_mean`: tracked `231.04444444444445`, fresh `255.7962962962963`
- `arms[0].frac_hitting_maxiter`: tracked `0.1`, fresh `0.18518518518518517`
- `arms[1].arm`: tracked `"fd"`, fresh `"ex"`
- `arms[1].solves`: tracked `90`, fresh `216`
- `arms[1].nit_mean`: tracked `159.32222222222222`, fresh `181.06944444444446`
- `arms[1].nfev_mean`: tracked `225.23333333333332`, fresh `255.7962962962963`
- `arms[1].frac_hitting_maxiter`: tracked `0.1`, fresh `0.18518518518518517`
- `arms[2].arm`: tracked `"an"`, fresh `"fd"`
- `arms[2].solves`: tracked `90`, fresh `216`
- `arms[2].nit_mean`: tracked `157.11111111111111`, fresh `180.65740740740742`

### project_stability / project_stability_partial77.json
- `agg.complete`: tracked `false`, fresh `true`
- `agg.requested`: tracked `126`, fresh `77`

