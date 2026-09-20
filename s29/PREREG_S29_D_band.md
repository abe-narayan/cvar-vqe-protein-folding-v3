# PREREG_S29_D_band -- WITHIN-REALISM-BAND ORDERING: THE ONE VERSION OF THE RECOGNITION QUESTION THE PERCEPTION-DISTORTION THEOREM LEAVES OPEN

Lane D (Adversary), Sprint 29. Written 2026-09-20 00:28 Pacific, **before the first number of
this experiment exists**. Contract `s29/S29_CONTRACT.md` rules 7, 10, 13, 16, 17, 18, 19 and 20.
Commissioned by the coordinator after lane L's S29-L12. Code `s29/s29_D_band.py`; results
`s29/results/s29_D_band*.json`; tests `tests/test_s29_D.py`. Addenda are appended, never edited,
once a result exists.

## 0. Labels
**EVERY NUMBER THIS EXPERIMENT PRODUCES IS ORACLE**: the ordering skill is measured against the
CA-RMSD to the native. Nothing here is deployable, nothing is tuned, no native quantity chooses
any parameter of any arm. The *consequence* of a positive would be a deployable two-stage
operator, and that operator would be a separate, separately pre-registered experiment.

## 1. The question, and why it is not already closed (rule 10)
Blau & Michaeli's theorem 3 (S29-L12) forbids a realism measure from preferring the
distortion-optimal answer **across** realism levels, and says nothing about ordering **inside**
one level. Every closure in the record -- S28-L35/L48 (20 of 31 scorers prefer production to a
0.25 A ORACLE structure), S28-L36 (the pool-member control: CAGEO prefers any real trace to a
contracted one), S21 L14/L17/L18 (in-pool selection exhausted), S12 (`nothing-ranks-within-the-pool`,
all-atom Amber reranking worth +0.004 A) -- measures the ACROSS-band question or measures
within-pool ordering **without conditioning on realism at all**. The new formulation is: condition
on a native-free realism statistic first, then ask for ordering skill inside the band. The
in-band metric of `in-band-is-the-only-ranking-metric` (S12) conditioned on the SCORE's own
top-24, which is a different conditioning (score band, not realism band) and is the thing
S28-L48's anti-correlation is measured on.
What would make this a re-run rather than a new question: if the realism statistic is
monotonically related to the scorer itself, "within band" degenerates to "within a narrow slice
of the scorer's own range" and the measurement is S12's in-band number again. **That is a real
risk and it is tested, not assumed** (section 5, the degeneracy check).

## 2. Findings engaged (rule 13)
Attacks finding 8 (recognition) in the one direction L's theorem leaves open; accepts 1 to 7 and
9 as binding; accepts 11 (common-mode error) as binding and irrelevant here (this is an ordering
question inside a fixed candidate set, not an aggregation question).

## 3. The instrument (reused, not re-implemented)
- Structures: for each of the 126 targets, (a) the **500 pool members** (`s27.run_pool.channels_for`
  via `s29.s29_D_cost_audit.load_target`), posed in lane A's production frame, and (b) the **nine
  ladder rungs** of the cost meter's cache (`s29/results/s29_D_ladder_structs/`), which carry
  PROD, the ORACLE near-native structures and the two matched controls.
- Scorers: the **31 S28-L48 scorers** through the meter's own adapters -- the 15 CA-level channels
  on CA clouds, and the 16 backbone channels, which need torsions and therefore are evaluated on
  the **projected chains** (the meter's `chain` basis; the projection is built once per structure
  and cached). Pool members carry their real torsions, so backbone scorers on pool members use
  the S27 cache (`s27/cache/<pdb>.npz`) exactly as S28-L36 did, and the cross-basis mismatch
  (projected ladder structures against real-torsion pool members) is stated beside every number,
  as S28-L48 stated it.
- ORACLE RMSD: `s12.instrument.ca_rmsd` to `nat_ca`, point cloud for CA-level scorers and built
  chain for the chain basis; both bases reported, the built chain deciding (contract rule 5).
- Statistics: `s24.stats_lib` (`compare`, `fmt`, `pinned_folds`); Spearman by `scipy.stats`;
  ties handled as in `ST.argmin_tied`'s spirit (a tied realism statistic never orders anything).

## 4. The realism statistic R, and the bands
**R is native-free and is NOT the scorer under test.** Three definitions, all pre-registered, run
as three arms (they are not a grid to pick from: each is reported, and the primary is R1):
- **R1 (primary), the pool-percentile of a REFERENCE realism channel.** R1(X) = the percentile of
  X's CAGEO value among the 500 real pool members' CAGEO values (`s27/cache`). CAGEO is the CA
  virtual-angle/torsion statistic: S28-L35/L36 showed it is the library's purest realism axis
  (it "prefers any real protein trace to a contracted one", h2h 0.962 vs a random signed
  combination) and that it carries no nativeness. It is therefore the right *conditioning*
  variable and the wrong *ranking* variable -- and when CAGEO itself is the scorer under test,
  that cell is reported as DEGENERATE and excluded from the headline.
- **R2, the geometric band.** R2(X) = the pair (mean virtual CA-CA bond, Rg) mapped to its
  percentile among the pool's, combined as the Euclidean distance to the pool's median in that
  two-dimensional percentile space. This is the realism axis the record has measured most often
  (contraction, S23 L1, `averaging-space-beats-the-objective`, S28-L36) and is scorer-independent.
- **R3, the leave-one-out consensus percentile.** R3(X) = the percentile of X's mean CA-RMSD to
  the 500 pool members among the members' own leave-one-out values: "is this structure as typical
  of the pool as a real member is". Scorer-independent; the same statistic `consensus-is-outlier-avoidance`
  priced (the native sits at its 82.8th percentile), which is precisely why it is a realism axis
  and not a nativeness axis.
**Bands**: the percentile range [0, 1] cut at the fixed edges {0.0, 0.2, 0.4, 0.6, 0.8, 1.0} (five
equal bands, fixed now, not chosen on any outcome). A structure's band is the band of its R value.
**Occupancy rule, fixed now**: a (target, band) cell contributes only if it holds **>= 20 distinct
structures with >= 10 distinct RMSD values**; cells below that are dropped and counted, and the
number of contributing cells is reported per scorer per band. A scorer that survives on fewer
than 60 of the 126 targets is reported as UNDERPOWERED and is not eligible for the headline.

## 5. The primary measurement, the controls, and the degeneracy check
For every scorer s, every band b, every target t with a contributing cell:
    rho_in(s, b, t) = Spearman( s(X), RMSD(X) ) over the structures X in cell (t, b).
Positive rho = the scorer orders by accuracy inside the band (lower score, lower RMSD).
Reported: the mean over targets per (s, b) with the **fold-clustered CI**, the pooled mean over
bands (weighting each target equally, bands averaged within target), the per-band occupancy, and:
- **The across-band comparison (the record's number).** rho_across(s, t) = Spearman over ALL
  structures of that target, which is the quantity S28-L48's ladder rho and S12's in-band metric
  are cousins of. The claim of this experiment is that rho_in > rho_across; the paired difference
  per target is reported with `ST.compare` and its fold CI.
- **The degeneracy check (mandatory, decides whether the arm means anything).** Per (s, b): the
  Spearman between s(X) and R(X) inside the band, and the interquartile range of s(X) inside the
  band as a fraction of its full-pool range. If |Spearman(s, R)| > 0.8 inside the band, or the
  band retains < 20% of the scorer's full range, the cell is DEGENERATE (the band is a slice of
  the scorer itself) and is excluded from the headline and reported separately.
- **The shuffle null.** Per (s, b), the same statistic with the RMSD values permuted WITHIN the
  cell (200 draws, seeded `SD.stable_rng(pdb, "s29D_band_null", salt="s29D")`): gives the null's
  mean and 95th percentile for rho_in at that cell's exact size, so a small-cell bias cannot
  masquerade as skill.
- **The max-over-scorers null (rule 17).** 31 scorers x 5 bands = 155 cells; the headline claim
  is priced against the distribution of the MAXIMUM mean rho_in over the 31 scorers under the
  per-target sign-flip null (`s27/s28_D_c2_chain_null.py`'s construction, transported), p_max
  reported. The band index is NOT a free parameter in the headline: the headline is the POOLED
  (band-averaged) number, and per-band numbers are reported as the decomposition.
- **The pool-only arm.** The ladder rungs are 9 structures against 500 pool members, so the
  headline is computed on POOL MEMBERS ONLY (a homogeneous set of real traces) and the version
  including the ladder rungs is reported beside it. This prevents the ORACLE rungs -- which are
  chosen against the native -- from carrying the correlation.

## 6. Falsifiers, registered now
- **F1 (the headline).** No scorer has pooled within-band rho_in with the fold-clustered CI
  excluding zero in the POSITIVE direction on the pool-only arm, after excluding DEGENERATE cells.
  If F1 fires, the recognition ceiling at this length is measured with controls and H0 gains its
  strongest evidence (a contribution, per the coordinator).
- **F2 (the mechanism).** Even if some scorer clears F1, the effect is not a recognition result
  unless rho_in > rho_across on the same targets with the fold CI of the paired difference
  excluding zero: otherwise the banding did nothing and the number is the record's own in-pool
  ordering re-measured.
- **F3 (the size that matters).** A positive pooled rho_in below 0.10 is reported as measured and
  is NOT called deployable: `in-band-ordering-is-per-target` puts the ordering needed for 2.0 A at
  0.638 and the achievable native-free proxies at 0.24 to 0.37, and `nothing-ranks-within-the-pool`
  prices perfect in-pool selection at 2.355 A against 3.324 A returned. The Angstrom consequence
  of any positive is computed by the ORACLE-within-band selection ceiling (below), never asserted.
- **F4 (the multiplicity).** A positive whose p_max over 31 scorers exceeds 0.05 is reported as
  the expected maximum of 31 nulls, as S28-L48's CONTACT@chain was.
**Registered prior (mine, stated before the run).** F1 does NOT fire for any scorer on the
pool-only arm: the record's in-pool ordering skill is at or below the noise on every channel
(S21 L14/L17/L18; S12 `nothing-ranks-within-the-pool`; S28-L48), and conditioning on realism
removes variance from the predictor without adding any signal about the native. I expect
rho_in ~ 0 to +0.05 pooled, with CAGEO's own cell DEGENERATE, and rho_in > rho_across simply
because rho_across is contaminated by the realism axis the theorem describes -- i.e. **F2 may
well pass while F1 fails, and that combination is not a result.** If anything surprises me it
will be DIS or DIS_MEAN in the middle band, where the distogram's own information is not
competing with a contraction signal.

## 7. The Angstrom consequence, computed not asserted
If any scorer clears F1 and F2, the same rows give the ORACLE diagnostic that prices it: the
mean RMSD of the structure the scorer picks INSIDE the band it lands production in, against
(a) production, (b) a random member of that band (the matched control: same band, no ordering),
(c) the band's ORACLE best (the ceiling of the two-stage operator). All three on the point cloud
and, for any cell that clears 0.7x MDE, on the built chain through `s12.instrument.project`.
That is a diagnostic; a deployable arm needs its own prereg, its own NaN-poison and its own
endpoint run.

## 8. Cost, jobs, multiplicity
One pass over 126 targets x 500 members x 31 scorers reuses the S27 channel cache for the pool
(no recomputation) and needs the 9 ladder rungs' scores, which the meter already computes. The
projected-chain basis for the 16 backbone scorers uses the ladder cache's chain rungs (built by
`s29_D_cost_audit.py build-cache --chain`, ~1 h) and the pool's real torsions from `s27/cache`.
Probe one target first and quote its peak RSS; then one governed job
(`--agent S29D --tag CPU --est-ram 1.0`), per-target checkpointed jsonl, resumable.
**Multiplicity declared: 0 endpoint comparisons** (nothing here touches the built-chain endpoint
of a deployable arm); 31 scorers x 5 bands x 3 realism definitions = 465 ORACLE diagnostic cells,
priced by the max-over-31 null on the primary (R1) and reported as a decomposition elsewhere;
the headline is ONE number per scorer (the pooled, band-averaged, pool-only rho_in).

## 9. What this experiment will not do
No deployable arm; no tuning of a band edge, a scorer weight or an occupancy threshold on any
outcome; no benchmark60; no AMBER; no second realism definition invented after seeing R1; no
re-cut of the bands after seeing a result (the five edges are fixed above); and no claim that a
positive within-band rho is an Angstrom until section 7's ceiling is computed.
