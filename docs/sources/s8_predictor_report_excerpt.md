# Source excerpt: `s8/predictor_report.log` (the evidence behind S8-13)

Written by Sprint 26, lane I (operational item 4), 2026-09-13. This file exists because the
log it quotes is the only on-disk source of a published number and the log itself is untracked
by decision (`.gitignore` treats `*.log` as disposable; `_archive/` is ignored as a whole).

## The file

| | |
|---|---|
| path on disk | `_archive/logs/s8/predictor_report.log` (moved from `s8/predictor_report.log` on 2026-09-04; see `_archive/README.txt`) |
| sha256 | `2647f640efb6d8f8e0398620647d503c366aa9a0209f17742d927366dd0f4740` |
| size | 13,588 bytes, 172 lines |
| writer | `s8/predictor.py` (console redirection; no code writes the path) |
| cited by | `docs/FINDINGS.md` section **S8-13** ("CORRECTION to S8-10: the 'transfer law' is the identity map where we operate"), line 3711 of the current file, which names it as `s8/predictor_report.log`; the same number is repeated in the corrections ledger at the top of `docs/FINDINGS.md` (the row on S8-10's transfer law) |
| tracked in git | no (only this excerpt is) |

## The number that has no other source

S8-13's central claim is the re-measured transfer law:

> **selected = 1.009 * ref + 0.011** (r = 0.981)

A search over every `.json`, `.md`, `.py`, `.txt` and `.log` file in the tree for the string
`1.009 * ref` finds exactly two files: `docs/FINDINGS.md` (the claim) and this log (its
evidence). No results JSON carries it. That line is log line 92, quoted with its block below.

S8-13 also quotes a second regression, `selected = 1.019 * ref - 0.001` (r = 0.984) "over
converged references". **That line is not in this log** and no other on-disk file contains it;
its source is not on this disk. Recorded here so the gap is known rather than discovered again.

## Verbatim excerpt 1: the transfer-law block (log lines 80 to 100)

```
--- sub-2 A and sub-1.5 A counts on the returned structure -----------------
  reference    <2.0  <1.5   (shipped selector: 27 and 20 of 126)
  O_gbest        54    40
  X_fit          35    23
  P_med          35    25
  P_fit          36    22
  Xw_fit         32    24
  P_fit2         37    22
  PfitB          33    25
  PfitS          33    22

--- TRANSFER LAW: does this file's references land on diffuse's line? -----
  over 2520 (target, reference) pairs:  selected = 1.009 * ref + 0.011   r = 0.981
  diffuse.py measured 0.940 * ref + 0.248 (r = 0.945) over 4,581 pairs
  the line crosses y = x at -1.23 A: BELOW that the arm loses

--- diagnostics -----------------------------------------------------------
  corr(AR log-lik of pool member, distogram score)           -0.306  (n=126)
  corr(AR log-lik of sample, distogram score)                -0.257  (n=126)
  CA-RMSD from the coordinate average to its projection      +1.050  (n=126)
  CA-CA step of the raw coordinate average (native 3.81)     +2.961  (n=126)
```

## Other S8-13 numbers that this log carries (log lines 1 to 5 and 102 to 119)

S8-13 also quotes, from the same run: X_fit converged at 3.217 (d -0.237 [-0.378, -0.098] vs
3.454, 79/47; -0.065 [-0.123, -0.006] vs S8-8's 3.282); P_fit 3.240; the three-filter ceiling
3.653 / 3.285 / 2.220 with the score worth -0.368 [-0.537, -0.204] and a perfect filter a
further -1.066 [-1.258, -0.885]. Those appear in `s8/consensus2_cache/` only as per-target
sweep records, not as the aggregated figures, so this log is their aggregated source too.

```

============================================================================================
S8-10  THE SINGLE-STRUCTURE PREDICTOR   n=126 targets
============================================================================================
  pool best 1.711   pool mean 4.453   shipped selector 3.454
```

```
--- THE PROJECTION, CONVERGED (8 sweeps / 36-point grid / 5 starts) -------
  n=126.  consensus2 reports 3.200 for `fit`; the build stage's 4/24/4 projection gave 3.285.
  arm        coarse  converged   gap A  pos-phi  d vs shipped             95% CI      W/L   d vs S8-8             95% CI
  P_fit       3.285      3.240   0.918   0.0648        -0.214 [-0.350,-0.083]   76/50       -0.042 [-0.102,+0.020]  60/66 
  X_fit       3.265      3.217   0.949   0.0629        -0.237 [-0.378,-0.098]   79/47       -0.065 [-0.123,-0.006]  74/52 
  Pws_fit     3.318      3.250   0.868   0.0705        -0.204 [-0.313,-0.100]   79/47       -0.032 [-0.107,+0.042]  67/59 
  P_med       3.282      3.282   0.000      nan        -0.172 [-0.316,-0.032]   75/51 
  the native instrument's positive-phi rate is 0.0558; `torfit` reached 0.1229

--- WHERE THE CEILING IS: same construction, three filters.  n=126 ---
  pool best 1.711   shipped selector 3.454
  filter     set mean  set best  CONSTRUCTED  pos-phi  d vs shipped             95% CI      W/L
  random        4.436     2.016        3.653   0.0445        +0.199 [-0.003,+0.407]   54/72 
  shipped       3.551     2.306        3.285   0.0517        -0.169 [-0.305,-0.037]   73/53 
  ORACLE        2.686     1.711        2.220   0.0498        -1.234 [-1.432,-1.041]  112/14 
  the shipped score is worth -0.368 [-0.537,-0.204] to the construction over no filter at all (84/42)
  a PERFECT filter would be worth a further -1.066 [-1.258,-0.885] (122/4) -- ORACLE
  and even then the construction stops at 2.220 A, against a pool best of 1.711 A
```

One discrepancy worth knowing: S8-13 states X_fit's positive-phi rate as 0.0551; the converged
X_fit row above reads 0.0629 (the 0.0551 may be the coarse projection's rate, which this table
does not print). Not resolved here; noted so nobody quotes either figure as sourced from this
log without checking.

## Verbatim excerpt 3: the ESM ablation (log lines 148 to 173), quoted by S8-13

```
--- THE ESM ABLATION: identical model, corpus, folds, schedule and RNG;
    the only difference is 32 input dimensions.  n=126 paired.
  estimator   ar (ESM)  ar_noesm  d = ESM - one-hot             95% CI      W/L
  best1 *        4.372     4.942             -0.570 [-0.891,-0.248]   81/45 
  best10 *       2.863     2.893             -0.030 [-0.188,+0.127]   60/66 
  best100 *      2.146     2.197             -0.051 [-0.153,+0.055]   68/58 
  medoid         3.846     4.007             -0.161 [-0.326,-0.013]   75/51 
  circmean       4.058     4.637             -0.578 [-0.894,-0.262]   83/43 
  llmax          3.942     4.240             -0.299 [-0.538,-0.075]   68/58 
  llmed          3.883     4.225             -0.342 [-0.572,-0.124]   73/53 
  llwc           3.883     4.348             -0.465 [-0.752,-0.189]   72/54 
  llfit          3.806     4.053             -0.247 [-0.414,-0.091]   76/50 
  dgmed          3.511     3.528             -0.017 [-0.145,+0.112]   64/62 
  dgfit          3.383     3.408             -0.025 [-0.127,+0.076]   67/59 
    * ORACLE rows: coverage of the ensemble, not a deployable estimate.
    negative d = ESM conditioning helps.

  validation NLL (nats/angle) on held-out PEPTIDES, per fold:
  fold    ar (ESM)  ar_noesm        d    best epoch (ESM / one-hot)
  0         0.7983    0.7715  +0.0268              7 /            9
  1         0.8455    0.8165  +0.0290              6 /           11
  2         0.9712    0.9315  +0.0398              4 /           11
  3         0.8049    0.7763  +0.0286              4 /           15
  4         0.8073    0.7986  +0.0087              9 /           12
  mean      0.8454    0.8189  +0.0266
```

## How to re-verify

```
sha256sum _archive/logs/s8/predictor_report.log
#   2647f640efb6d8f8e0398620647d503c366aa9a0209f17742d927366dd0f4740
grep -n "selected = 1.009" _archive/logs/s8/predictor_report.log docs/FINDINGS.md
```

`python s26/examine.py` runs this check as claim `s8_13_transfer_law` and reports ABSENT (with
this excerpt as the fallback) if the archive directory is ever lost.
