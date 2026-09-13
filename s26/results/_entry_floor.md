`s26/ph_cis.py floor`, `s26/results/ph_cis_floor.json` (complete 126/126), job
`s26/jobs_done/ph_cis_floor.json` (exit 0, 10 s, peak RSS 0.098 GB, held 345 s at the job cap).
Pre-registered in `s26/PREREG_cis.md` part 2 and addendum 1. ORACLE DIAGNOSTIC: every number
reads the native. Bases, named: `floor_ca` is the CA-RMSD between the model-1 native CA trace
and the chain rebuilt from the native's OWN phi/psi at ideal trans geometry
(`core.geometry.build_backbone`), after superposition; `rebuild_bb` is the same on N/CA/C
(`core.data.Peptide.rebuild`, verified from `core/data.py:410-419`); `chain_cost` is the
production `rmsd_arm - rmsd_avg` (built chain minus point cloud, a difference across bases,
`bench_results/cache/1fc9f2dcf489e2fb`).

    floor_ca   (CA, ideal trans on the native's own torsions)   mean 0.347 A (SE 0.029), median 0.272, p90 0.752, max 1.474 (1ID6)
    rebuild_bb (N/CA/C, `Peptide.rebuild`)                         mean 0.340 A (SE 0.028), median 0.248, max 1.413
    chain_cost (production rmsd_arm - rmsd_avg)                    mean 0.166 A (SE 0.018), median 0.098
    targets with floor_ca > 0.5 A                                  32 / 126;  floor_ca > chain_cost on 86 / 126
    Spearman(floor_ca, max omega deviation)                        +0.828;  with mean omega deviation +0.833
    Spearman(chain_cost, max omega deviation)                      -0.036;  Spearman(chain_cost, floor_ca) +0.083
    cis targets                                                    0 (L22); the cis-vs-non-cis contrast is empty and is not printed
    9UV5 (the one bond beyond 30 deg)                              floor_ca 0.663, rebuild_bb 0.708, chain_cost 0.372

`ST.fmt`, verbatim (the two quantities are both per-target Angstroms but are different
objects; "BETTER" below means only that the production chain cost is smaller than the
own-torsion rebuild floor, not that anything improved):

    production chain cost (rmsd_arm - rmsd_avg) MINUS the ORACLE CA floor of ideal trans geometry on the native's own torsions
      a 0.1664 (med 0.0977)   b 0.3468 (med 0.2722)   n=126
      effect -0.1804   median -0.1261   SE 0.0333   MDE 0.0933   effect/MDE -1.93
      iid  CI95 [-0.2433, -0.1176]
      fold CI95 [-0.2266, -0.1022]   folds same sign 5/5   per-fold 0:-0.025 1:-0.224 2:-0.230 3:-0.223 4:-0.202
      86W/40L/0T   worst degradation +0.5161 (1RSW)   p90 +0.2720   power 1.00  Type-M 1.00
      concentration: drop-top10 -0.1100 vs uniform-effect null p10/p50/p90 -0.1523/-0.1116/-0.0702 -> pctile 0.518
      VERDICT: BETTER

Reading. (1) The constant omega does carry a representation cost on this instrument even
without a single cis bond: the ideal-trans rebuild of the native's own torsions misses the
native by 0.35 A on average and by more than 0.5 A on a quarter of the targets, and that miss is
the omega non-planarity to rho 0.83 (L22: 0.5% of bonds beyond 20 deg, mean deviation 1.9 deg;
small deviations accumulate along the chain). (2) It is NOT what the production projection
pays: the chain cost (0.166 A) does not correlate with the floor (rho +0.08) nor with the omega
deviation (rho -0.04). The projection's cost is a displacement effect of the operator (S16 L27:
the projection moves 0.813 A with a negative cosine), not a representation effect. (3)
`floor_ca` is an UPPER bound on the manifold floor, because the projection fits phi/psi to the
trace rather than rebuilding from the native's torsions; `s26/PREREG_cis.md` addendum 2
registers `floor2` (the native projected through the production projection) as the tight
number, 10 min CPU, to run after the reject jobs. Power: n = 126, SE 0.033, MDE 0.093 on the
contrast; a Spearman at n = 126 has SE about 0.09, so the two nulls (+0.08, -0.04) exclude
|rho| above about 0.25. Nothing here is a proposal; the two-bond-length projection stays worth
0.000 A on this instrument (L22), and the design note in `s26/agentPH_FINDINGS.md` 1.4 stands.
